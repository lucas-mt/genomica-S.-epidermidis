"""
#na linha de comando
python3 strains_with_virus <arq>.fasta
"""

import argparse
from pathlib import Path
from collections import Counter


def contar_accessions(caminho_fasta):

    caminho_fasta = Path(caminho_fasta)

    if not caminho_fasta.exists():
        raise FileNotFoundError(
            f"Arquivo não encontrado: {caminho_fasta}"
        )

    contagem = Counter()
    total_cabecalhos = 0

    print("Iniciando leitura do FASTA...")

    with open(caminho_fasta, "r", encoding="utf-8") as fasta:

        for linha in fasta:

            if linha.startswith(">"):

                total_cabecalhos += 1

                partes = linha.split("accn|", 1)

                if len(partes) < 2:
                    continue

                accession = partes[1].split("|", 1)[0].strip()

                if accession:
                    contagem[accession] += 1

                if total_cabecalhos % 100000 == 0:
                    print(
                        f"Cabeçalhos processados: {total_cabecalhos}"
                    )

    print("\nLeitura concluída!")
    print(f"Total de cabeçalhos: {total_cabecalhos}")
    print(f"Accessions diferentes: {len(contagem)}")

    return contagem


def salvar_resultado(contagem, caminho_saida):

    print("\nSalvando resultado...")

    with open(caminho_saida, "w", encoding="utf-8") as arquivo:

        arquivo.write("accession\tquantidade\n")

        for accession, quantidade in sorted(contagem.items(), key=lambda item: item[1], reverse=True):

            arquivo.write(f"{accession}\t{quantidade}\n")

    print(f"Arquivo salvo em: {caminho_saida}")


def main():

    parser = argparse.ArgumentParser(
        description="Conta quantas vezes cada accession aparece em um FASTA."
    )

    parser.add_argument(
        "fasta",
        help="Caminho do arquivo FASTA."
    )

    parser.add_argument(
        "-o", "--output",
        default="contagem_accessions.tab",
        help="Nome do arquivo TAB de saída."
    )

    args = parser.parse_args()

    contagem = contar_accessions(args.fasta)

    salvar_resultado(contagem, args.output)


if __name__ == "__main__":
    main()
