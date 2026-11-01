"""
#na linha de comando
python3 verificar_bases.py original.fasta cromossomos.fasta plasmideos.fasta
"""

import argparse
from pathlib import Path
from collections import defaultdict
from openpyxl import Workbook


def ler_fasta(caminho):
    """
    Lê um arquivo FASTA e retorna cada cabeçalho
    acompanhado da quantidade de bases da sequência.
    """

    cabecalho = None
    bases = 0

    with open(caminho, "r", encoding="utf-8", errors="replace") as arquivo:

        for linha in arquivo:

            linha = linha.strip()

            if not linha:
                continue

            if linha.startswith(">"):

                if cabecalho is not None:
                    yield cabecalho, bases

                cabecalho = linha
                bases = 0

            else:
                bases += len(linha)

        if cabecalho is not None:
            yield cabecalho, bases


def carregar_original(caminho):
    """
    Armazena os cabeçalhos e suas respectivas quantidades de bases.
    """

    registros = defaultdict(list)

    for cabecalho, bases in ler_fasta(caminho):
        registros[cabecalho].append(bases)

    return registros


def verificar_fasta(original, destino, nome_aba):
    """
    Compara a quantidade de bases dos registros do FASTA
    de destino com os registros correspondentes do original.
    """

    resultados = []

    ocorrencias = defaultdict(int)

    for cabecalho_destino, bases_destino in ler_fasta(destino):

        ocorrencias[cabecalho_destino] += 1
        numero = ocorrencias[cabecalho_destino]

        bases_originais = original.get(cabecalho_destino, [])

        if numero > len(bases_originais):

            cabecalho_original = ""
            bases_original = None
            diferenca = None
            status = "CABECALHO_NAO_ENCONTRADO"

        else:

            cabecalho_original = cabecalho_destino
            bases_original = bases_originais[numero - 1]

            diferenca = bases_destino - bases_original

            if diferenca == 0:
                status = "OK"
            else:
                status = "DIFERENCA_DE_BASES"

        resultados.append([
            len(resultados) + 1,
            cabecalho_original,
            bases_original,
            cabecalho_destino,
            bases_destino,
            diferenca,
            status
        ])

    return resultados


def criar_planilha(resultados_cromossomos, resultados_plasmideos, saida):

    workbook = Workbook()

    abas = {
        "cromossomos": resultados_cromossomos,
        "plasmideos": resultados_plasmideos
    }

    for nome_aba, resultados in abas.items():

        planilha = workbook.create_sheet(nome_aba)

        cabecalhos = [
            "registro",
            "cabecalho_original",
            "bases_original",
            "cabecalho_destino",
            "bases_destino",
            "diferenca_bases",
            "status"
        ]

        planilha.append(cabecalhos)

        for registro in resultados:
            planilha.append(registro)

        planilha.freeze_panes = "A2"
        planilha.auto_filter.ref = planilha.dimensions

        larguras = {
            "A": 12,
            "B": 75,
            "C": 18,
            "D": 75,
            "E": 18,
            "F": 18,
            "G": 30
        }

        for coluna, largura in larguras.items():
            planilha.column_dimensions[coluna].width = largura

    # Remove a aba padrão criada pelo openpyxl
    if "Sheet" in workbook.sheetnames:
        del workbook["Sheet"]

    workbook.save(saida)


def main():

    parser = argparse.ArgumentParser(
        description="Compara a quantidade de bases entre arquivos FASTA."
    )

    parser.add_argument(
        "original",
        help="Caminho do FASTA original."
    )

    parser.add_argument(
        "cromossomos",
        help="Caminho do FASTA contendo os cromossomos."
    )

    parser.add_argument(
        "plasmideos",
        help="Caminho do FASTA contendo os plasmídeos."
    )

    parser.add_argument(
        "-o",
        "--saida",
        default="verificacao_bases.xlsx",
        help="Nome do arquivo Excel de saída."
    )

    args = parser.parse_args()

    arquivos = [
        args.original,
        args.cromossomos,
        args.plasmideos
    ]

    for arquivo in arquivos:
        if not Path(arquivo).is_file():
            parser.error(f"Arquivo não encontrado: {arquivo}")

    print("\nLendo FASTA original...")

    original = carregar_original(args.original)

    print(f"Registros encontrados: {sum(len(v) for v in original.values())}")

    print("\nVerificando cromossomos...")

    resultados_cromossomos = verificar_fasta(
        original,
        args.cromossomos,
        "cromossomos"
    )

    print(f"Registros verificados: {len(resultados_cromossomos)}")

    print("\nVerificando plasmídeos...")

    resultados_plasmideos = verificar_fasta(
        original,
        args.plasmideos,
        "plasmideos"
    )

    print(f"Registros verificados: {len(resultados_plasmideos)}")

    print("\nGerando planilha de verificação...")

    criar_planilha(
        resultados_cromossomos,
        resultados_plasmideos,
        args.saida
    )

    print("\nVerificação concluída!")
    print(f"Arquivo gerado: {args.saida}")


if __name__ == "__main__":
    main()
