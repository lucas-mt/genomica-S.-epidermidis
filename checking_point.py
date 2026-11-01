"""
#na linha de comando
python3 checking_point <chrom/plasm>.ods <cromossomos/plasmideos>.fasta
"""

import argparse
import re
import sys
from pathlib import Path

import pandas as pd


def normalizar_accession(valor):
    """Normaliza accessions para comparação."""
    if pd.isna(valor):
        return ""

    return str(valor).strip().upper()


def extrair_accession(cabecalho):
    """
    Extrai o accession do cabeçalho FASTA.

    Formato principal:
    >accn|CP118842| Staphylococcus epidermidis ...

    Também tenta identificar accessions quando o prefixo
    accn| não está presente.
    """
    # Primeiro, procura o formato padrão.
    correspondencia = re.match(
        r"^>accn\|([^|\s]+)\|",
        cabecalho.strip(),
        flags=re.IGNORECASE
    )

    if correspondencia:
        return normalizar_accession(correspondencia.group(1))

    # Alternativamente, procura accessions GenBank
    # em outros formatos de cabeçalho.
    padrao = (
        r"(?<![A-Za-z0-9])"
        r"([A-Z]{1,4}_?\d{5,}(?:\.\d+)?)"
        r"(?![A-Za-z0-9])"
    )

    correspondencia = re.search(
        padrao,
        cabecalho,
        flags=re.IGNORECASE
    )

    if correspondencia:
        return normalizar_accession(correspondencia.group(1))

    return ""


def ler_planilha(caminho):
    """Lê os accessions da aba main."""
    caminho = Path(caminho)

    if not caminho.exists():
        raise FileNotFoundError(
            f"Planilha não encontrada: {caminho}"
        )

    df = pd.read_excel(
        caminho,
        sheet_name="main",
        dtype=str
    )

    # Localiza a coluna accession sem depender de
    # maiúsculas, minúsculas ou espaços extras.
    colunas = {
        str(coluna).strip().lower(): coluna
        for coluna in df.columns
    }

    if "accession" not in colunas:
        raise ValueError(
            "A coluna 'accession' não foi encontrada "
            f"na aba 'main'. Colunas disponíveis: "
            f"{list(df.columns)}"
        )

    coluna_accession = colunas["accession"]

    # Mantém o número da linha original da planilha.
    df = df.copy()
    df["_linha_planilha"] = range(2, len(df) + 2)
    df["_accession_normalizado"] = (
        df[coluna_accession].apply(normalizar_accession)
    )

    return df, coluna_accession


def ler_cabecalhos_fasta(caminho):
    """
    Lê somente os cabeçalhos FASTA.
    Não carrega as sequências na memória.
    """
    caminho = Path(caminho)

    if not caminho.exists():
        raise FileNotFoundError(
            f"Arquivo FASTA não encontrado: {caminho}"
        )

    registros = []

    with caminho.open(
        "r",
        encoding="utf-8",
        errors="replace"
    ) as arquivo:

        for numero_linha, linha in enumerate(arquivo, 1):
            linha = linha.rstrip("\r\n")

            if linha.startswith(">"):
                accession = extrair_accession(linha)

                registros.append({
                    "fasta_registro": len(registros) + 1,
                    "linha_fasta": numero_linha,
                    "fasta_accession": accession,
                    "fasta_header": linha,
                })

    return pd.DataFrame(
        registros,
        columns=[
            "fasta_registro",
            "linha_fasta",
            "fasta_accession",
            "fasta_header",
        ]
    )


def comparar(planilha, fasta, saida):
    """Compara os accessions e gera uma nova planilha."""

    df, coluna_accession = ler_planilha(planilha)
    df_fasta = ler_cabecalhos_fasta(fasta)

    # Considera apenas accessions preenchidos na planilha.
    df_validos = df[
        df["_accession_normalizado"] != ""
    ].copy()

    # Accessions únicos para identificar os registros
    # que não possuem correspondência.
    accessions_planilha = set(
        df_validos["_accession_normalizado"]
    )

    accessions_fasta = set(
        df_fasta.loc[
            df_fasta["fasta_accession"] != "",
            "fasta_accession"
        ]
    )

    # --------------------------------------------------
    # 1. Correspondências
    # --------------------------------------------------

    # Um accession pode aparecer mais de uma vez na
    # planilha ou no FASTA. Todas as ocorrências serão
    # preservadas para facilitar a conferência.
    df_comparacao = df_validos.merge(
        df_fasta[
            df_fasta["fasta_accession"] != ""
        ],
        left_on="_accession_normalizado",
        right_on="fasta_accession",
        how="inner"
    )

    correspondencias = pd.DataFrame({
        "linha_planilha": df_comparacao["_linha_planilha"],
        "accession_planilha": df_comparacao[coluna_accession],
        "accession_fasta": df_comparacao["fasta_accession"],
        "linha_fasta": df_comparacao["linha_fasta"],
        "registro_fasta": df_comparacao["fasta_registro"],
        "fasta_header": df_comparacao["fasta_header"],
    })

    # Identifica correspondências repetidas.
    contagem_planilha = (
        df_validos["_accession_normalizado"]
        .value_counts()
        .to_dict()
    )

    contagem_fasta = (
        df_fasta.loc[
            df_fasta["fasta_accession"] != "",
            "fasta_accession"
        ]
        .value_counts()
        .to_dict()
    )

    correspondencias["duplicado_planilha"] = (
        correspondencias["accession_planilha"]
        .map(normalizar_accession)
        .map(contagem_planilha)
        .gt(1)
    )

    correspondencias["duplicado_fasta"] = (
        correspondencias["accession_fasta"]
        .map(contagem_fasta)
        .gt(1)
    )

    # --------------------------------------------------
    # 2. Accessions presentes apenas na planilha
    # --------------------------------------------------

    apenas_planilha = df[
        (df["_accession_normalizado"] == "") |
        (~df["_accession_normalizado"].isin(accessions_fasta))
    ].copy()

    apenas_planilha = apenas_planilha.drop(
        columns=["_accession_normalizado"],
        errors="ignore"
    )

    apenas_planilha = apenas_planilha.rename(
        columns={
            coluna_accession: "accession",
            "_linha_planilha": "linha_planilha",
        }
    )

    # --------------------------------------------------
    # 3. Registros presentes apenas no FASTA
    # --------------------------------------------------

    apenas_fasta = df_fasta[
        (df_fasta["fasta_accession"] == "") |
        (~df_fasta["fasta_accession"].isin(accessions_planilha))
    ].copy()

    apenas_fasta["motivo"] = apenas_fasta[
        "fasta_accession"
    ].apply(
        lambda x: (
            "Cabeçalho sem accession identificável"
            if not x
            else "Accession não encontrado na planilha"
        )
    )

    # --------------------------------------------------
    # 4. Exportação
    # --------------------------------------------------

    caminho_saida = Path(saida)

    if caminho_saida.resolve() in {
        Path(planilha).resolve(),
        Path(fasta).resolve(),
    }:
        raise ValueError(
            "O arquivo de saída deve ser diferente "
            "dos arquivos de entrada."
        )

    with pd.ExcelWriter(
        caminho_saida,
        engine="openpyxl"
    ) as writer:

        correspondencias.to_excel(
            writer,
            sheet_name="correspondencias",
            index=False
        )

        apenas_planilha.to_excel(
            writer,
            sheet_name="so_planilha",
            index=False
        )

        apenas_fasta.to_excel(
            writer,
            sheet_name="so_fasta",
            index=False
        )

    print("\nComparação concluída!")
    print(f"Linhas da planilha: {len(df)}")
    print(f"Cabeçalhos FASTA: {len(df_fasta)}")
    print(f"Correspondências: {len(correspondencias)}")
    print(f"Registros apenas na planilha: {len(apenas_planilha)}")
    print(f"Registros apenas no FASTA: {len(apenas_fasta)}")
    print(f"\nArquivo gerado: {caminho_saida.resolve()}")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Compara accessions da aba 'main' de uma "
            "planilha com os cabeçalhos de um arquivo FASTA."
        )
    )

    parser.add_argument(
        "planilha",
        help="Caminho da planilha .ods ou .xlsx"
    )

    parser.add_argument(
        "fasta",
        help="Caminho do arquivo FASTA"
    )

    parser.add_argument(
        "-o",
        "--saida",
        default="comparacao_accessions.xlsx",
        help=(
            "Nome do arquivo Excel de saída "
            "(padrão: comparacao_accessions.xlsx)"
        )
    )

    args = parser.parse_args()

    try:
        comparar(
            args.planilha,
            args.fasta,
            args.saida
        )

    except Exception as erro:
        print(f"ERRO: {erro}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()

