import pandas as pd
import argparse
from pathlib import Path


def normalizar(valor):
    """Normaliza textos para facilitar a correspondência."""
    if pd.isna(valor):
        return None

    return str(valor).strip()


def processar_planilha(caminho):
    caminho = Path(caminho)

    if not caminho.exists():
        raise FileNotFoundError(
            f"Arquivo não encontrado: {caminho}"
        )

    if caminho.suffix.lower() != ".ods":
        raise ValueError(
            "O arquivo precisa estar no formato .ods."
        )

    # Lê todas as abas como texto.
    abas = pd.read_excel(
        caminho,
        sheet_name=None,
        engine="odf",
        dtype=str
    )

    # Verifica as abas necessárias.
    for nome in ["cepas com vírus", "original", "main"]:
        if nome not in abas:
            raise ValueError(
                f"A aba '{nome}' não foi encontrada."
            )

    virus = abas["cepas com vírus"].copy()
    original = abas["original"].copy()
    main = abas["main"].copy()

    # Remove espaços extras dos nomes das colunas.
    virus.columns = virus.columns.str.strip()
    original.columns = original.columns.str.strip()
    main.columns = main.columns.str.strip()

    # Verifica as colunas necessárias.
    colunas_virus = [
        "Genome ID",
        "Quantidade de vezes"
    ]

    colunas_original = [
        "Genome ID",
        "Strain"
    ]

    colunas_main = [
        "Strain",
        "accession"
    ]

    for coluna in colunas_virus:
        if coluna not in virus.columns:
            raise ValueError(
                f"Coluna '{coluna}' ausente na aba "
                "'cepas com vírus'."
            )

    for coluna in colunas_original:
        if coluna not in original.columns:
            raise ValueError(
                f"Coluna '{coluna}' ausente na aba 'original'."
            )

    for coluna in colunas_main:
        if coluna not in main.columns:
            raise ValueError(
                f"Coluna '{coluna}' ausente na aba 'main'."
            )

    # Normaliza os Genome IDs e os nomes das cepas.
    virus["_id_busca"] = virus["Genome ID"].map(normalizar)

    original["_id_busca"] = original["Genome ID"].map(
        normalizar
    )

    original["_Strain_busca"] = original["Strain"].map(
        normalizar
    )

    main["_Strain_busca"] = main["Strain"].map(
        normalizar
    )

    # Recupera o Strain correspondente a cada Genome ID.
    dados_Strain = original[
        ["_id_busca", "_Strain_busca"]
    ].dropna(subset=["_id_busca"])

    # Evita duplicar os registros da tabela de consulta.
    dados_Strain = dados_Strain.drop_duplicates(
        subset="_id_busca",
        keep="first"
    )

    mapa_Strain = dados_Strain.set
