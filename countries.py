import pandas as pd
import argparse
from pathlib import Path


def analisar_paises(caminho_planilha):

    caminho = Path(caminho_planilha)

    if not caminho.exists():
        raise FileNotFoundError(
            f"Arquivo não encontrado: {caminho}"
        )

    # Lê a aba main
    df = pd.read_excel(
        caminho,
        sheet_name="main",
        engine="odf"
    )

    coluna = "geo_loc_name"

    if coluna not in df.columns:
        raise ValueError(
            f"A coluna '{coluna}' não foi encontrada na aba main."
        )

    # Verifica se existem genomas
    if df.empty:
        raise ValueError("A planilha não contém genomas.")

    def normalizar(valor):

        if pd.isna(valor):
            return None

        valor = str(valor).strip()

        # Identifica células vazias ou contendo apenas "!"
        if valor == "" or valor == "!":
            return None

        # Padroniza maiúsculas e minúsculas
        return valor.casefold()


    paises = df[coluna].apply(normalizar)

    # Identifica informações ausentes
    paises = paises.fillna(
        "Country information unavailable"
    )

    # Conta os genomas por país
    contagem = paises.value_counts()

    # Calcula as porcentagens
    porcentagens = (contagem / len(df)) * 100

    # Cria a tabela final
    resultado = pd.DataFrame({
        "Country": porcentagens.index,
        "Percentage of genomes": porcentagens.values
    })

    # Padroniza a apresentação dos nomes
    resultado["Country"] = resultado["Country"].apply(
        lambda x: x.title()
        if x != "Country information unavailable"
        else x
    )

    # Arredonda as porcentagens
    resultado["Percentage of genomes"] = (
        resultado["Percentage of genomes"].round(2)
    )

    # Separa os países com informação dos sem informação
    sem_informacao = resultado[
        resultado["Country"] == "Country information unavailable"
    ]

    com_informacao = resultado[
        resultado["Country"] != "Country information unavailable"
    ]

    # Ordena os países do maior para o menor percentual
    com_informacao = com_informacao.sort_values(
        by="Percentage of genomes",
        ascending=False
    )

    # Coloca os países sem informação sempre no final
    resultado = pd.concat(
        [com_informacao, sem_informacao],
        ignore_index=True
    )

    # Define o nome do arquivo de saída
    saida = caminho.with_name(
        caminho.stem + "_country_distribution.tab"
    )

    # Exporta como arquivo tabular
    resultado.to_csv(
        saida,
        index=False,
        sep="\t",
        decimal=".",
        encoding="utf-8"
    )

    # Exibe a tabela no terminal
    print("\nDistribution of genomes by country:\n")
    print(resultado.to_string(index=False))

    print(f"\nResultado salvo em: {saida}")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Calculate genome distribution by country."
    )

    parser.add_argument(
        "planilha",
        help="Path to the .ods spreadsheet"
    )

    args = parser.parse_args()

    analisar_paises(args.planilha)

