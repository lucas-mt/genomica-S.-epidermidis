import argparse
import pandas as pd
import time

from Bio import Entrez, SeqIO
from io import StringIO
from datetime import datetime
from pathlib import Path
from tqdm import tqdm


# ==========================================================
# CONFIGURAÇÕES
# ==========================================================

Entrez.email = "lucas.lira@ufv.br"


# ==========================================================
# FUNÇÕES DE CONSULTA
# ==========================================================

def lista_accessions(planilha, tipo="odf"):

    db = pd.read_excel(planilha, engine=tipo)

    if "GenBank Accessions" not in db.columns:
        raise ValueError(
            'A coluna "GenBank Accessions" não foi encontrada.'
        )

    return list(db["GenBank Accessions"])


def IdList(accession):

    handle = Entrez.esearch(
        db="nuccore",
        term=f"{accession}[Accession]",
        retmax=1
    )

    resultado = Entrez.read(handle)
    handle.close()

    ids = resultado["IdList"]

    if not ids:
        raise ValueError("Accession não encontrado no NCBI.")

    return ids[0]


def infos_gerais(ID):

    handle = Entrez.efetch(
        db="nuccore",
        id=ID,
        rettype="gb",
        retmode="text"
    )

    texto = handle.read()
    handle.close()

    return texto


def trans_info(info, tipo="genbank"):

    registro = SeqIO.read(StringIO(info), tipo)

    return registro


def resolve_none(funcao, info):

    dado = funcao(info)

    if dado:
        return dado[0]

    return None


def coleta_ano(data):

    if not data:
        return None

    for i in data:

        if len(i) == 4 and i.isdigit():
            return i

    return None


def dados(accession, registro):

    resultado = {}

    for feature in registro.features:

        if feature.type == "source":

            qualificadores = feature.qualifiers

            collection_date = qualificadores.get(
                "collection_date"
            )

            geo_loc_name = qualificadores.get(
                "geo_loc_name"
            )

            referencias = registro.annotations.get(
                "references", []
            )

            referencia = (
                referencias[0]
                if referencias else None
            )

            data_registro = registro.annotations.get("date")

            resultado = {

                "organism": resolve_none(
                    qualificadores.get, "organism"
                ),

                "accession": accession,

                "strain": resolve_none(
                    qualificadores.get, "strain"
                ),

                "isolation_source": resolve_none(
                    qualificadores.get, "isolation_source"
                ),

                "geo_loc_name": (
                    geo_loc_name[0].split(":")[0].strip()
                    if geo_loc_name else None
                ),

                "collection_date": (
                    coleta_ano(
                        collection_date[0].split("-")
                    )
                    if collection_date else None
                ),

                "host": resolve_none(
                    qualificadores.get, "host"
                ),

                "registration_date": (
                    datetime.strptime(
                        data_registro, "%d-%b-%Y"
                    ).year
                    if data_registro else None
                ),

                "title": (
                    referencia.title
                    if referencia else None
                ),

                "journal": (
                    referencia.journal
                    if referencia else None
                )
            }

            break

    if not resultado:
        raise ValueError(
            "Registro não contém uma feature source."
        )

    return resultado


def dados_preenchidos(accession):

    ID = IdList(accession)

    infos = trans_info(
        infos_gerais(ID)
    )

    return dados(accession, infos)


# ==========================================================
# EXECUÇÃO PRINCIPAL
# ==========================================================

def main():

    parser = argparse.ArgumentParser(
        description=(
            "Consulta metadados de accessions no NCBI "
            "a partir de uma planilha."
        )
    )

    parser.add_argument(
        "planilha",
        help="Caminho da planilha (.ods, .xlsx ou .xls)."
    )

    parser.add_argument(
        "-o",
        "--saida",
        help="Nome do arquivo de saída (.tab)."
    )

    args = parser.parse_args()

    planilha = Path(args.planilha).expanduser()

    if not planilha.exists():
        parser.error(
            f"Arquivo não encontrado: {planilha}"
        )

    motores = {
        ".ods": "odf",
        ".xlsx": "openpyxl",
        ".xls": "xlrd"
    }

    formato = motores.get(planilha.suffix.lower())

    if formato is None:
        parser.error(
            "Formato não suportado. "
            "Utilize .ods, .xlsx ou .xls."
        )

    # Define o arquivo de saída.

    if args.saida:

        arquivo_saida = Path(args.saida).expanduser()

    else:

        arquivo_saida = planilha.with_name(
            planilha.stem + "_corrigido.tab"
        )

    # ======================================================
    # LEITURA DOS ACCESSIONS
    # ======================================================

    try:

        accessions = lista_accessions(
            planilha,
            formato
        )

    except Exception as erro:

        parser.error(
            f"Erro ao ler a planilha: {erro}"
        )

    # Expande células com múltiplos accessions.

    lista = []

    for accession in accessions:

        if pd.isna(accession):
            continue

        for acc in str(accession).split(","):

            acc = acc.strip()

            if acc:
                lista.append(acc)

    # Remove duplicatas preservando a ordem.

    lista = list(dict.fromkeys(lista))

    total = len(lista)

    if total == 0:

        print("Nenhum accession encontrado na planilha.")
        return

    print("\n" + "=" * 60)
    print("CONSULTA DE METADADOS NO NCBI")
    print("=" * 60)

    print(f"\nPlanilha: {planilha.name}")
    print(f"Total de accessions: {total}")
    print(f"Arquivo de saída: {arquivo_saida.name}")

    print("\nIniciando processamento...\n")

    # ======================================================
    # CONSULTAS COM BARRA DE PROGRESSO
    # ======================================================

    infos = []
    falhas = []

    inicio = time.time()

    for i, accession in enumerate(
        tqdm(
            lista,
            desc="Consultando NCBI",
            unit="accession",
            dynamic_ncols=True
        ),
        start=1
    ):

        try:

            resultado = dados_preenchidos(accession)

            infos.append(resultado)

            status = "OK"

        except Exception as erro:

            falhas.append({
                "accession": accession,
                "erro": str(erro)
            })

            status = "ERRO"

        # Atualiza as informações exibidas na barra.

        tqdm.write(
            f"[{i}/{total}] {accession}: {status}"
        )

        # Respeita o intervalo entre consultas.

        if i < total:
            time.sleep(0.34)

    # ======================================================
    # SALVAMENTO DOS RESULTADOS
    # ======================================================

    print("\nSalvando resultados...")

    colunas = [
        "organism",
        "accession",
        "strain",
        "isolation_source",
        "geo_loc_name",
        "collection_date",
        "host",
        "registration_date",
        "title",
        "journal"
    ]

    df = pd.DataFrame(
        infos,
        columns=colunas
    )

    df.to_csv(
        arquivo_saida,
        sep="\t",
        index=False,
        encoding="utf-8"
    )

    # Salva os accessions que apresentaram erros.

    if falhas:

        arquivo_falhas = arquivo_saida.with_name(
            arquivo_saida.stem + "_falhas.txt"
        )

        pd.DataFrame(falhas).to_csv(
            arquivo_falhas,
            sep="\t",
            index=False,
            encoding="utf-8"
        )

        print(
            f"\nArquivo de falhas: {arquivo_falhas}"
        )

    # ======================================================
    # RELATÓRIO FINAL
    # ======================================================

    tempo_total = time.time() - inicio

    minutos = int(tempo_total // 60)
    segundos = int(tempo_total % 60)

    print("\n" + "=" * 60)
    print("PROCESSAMENTO CONCLUÍDO!")
    print("=" * 60)

    print(f"Total de accessions: {total}")
    print(f"Consultas bem-sucedidas: {len(infos)}")
    print(f"Consultas com falha: {len(falhas)}")

    print(
        f"Tempo total: {minutos} min {segundos} s"
    )

    print(
        f"\nArquivo salvo em: {arquivo_saida.resolve()}"
    )

    print("=" * 60)


# ==========================================================
# INICIALIZAÇÃO
# ==========================================================

if __name__ == "__main__":

    main()
