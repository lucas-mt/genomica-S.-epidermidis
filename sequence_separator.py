import pandas as pd
import sys
import os

# Verificar se o nome da planilha foi informado
if len(sys.argv) != 2:
    print("Uso: python3 separar_cepas.py <planilha.ods>")
    sys.exit(1)

# Obter o nome da planilha pelo terminal
arquivo_entrada = sys.argv[1]

# Verificar se o arquivo existe
if not os.path.isfile(arquivo_entrada):
    print(f"Erro: o arquivo '{arquivo_entrada}' não foi encontrado.")
    sys.exit(1)

# Criar os nomes dos arquivos de saída automaticamente
nome, extensao = os.path.splitext(arquivo_entrada)

arquivo_cromossomos = f"{nome}_cromossomos.tab"
arquivo_plasmidios = f"{nome}_plasmidios.tab"

# Ler a aba "tudo"
df = pd.read_excel(
    arquivo_entrada,
    sheet_name="tudo",
    engine="odf"
)

# Verificar se a coluna necessária existe
if "NCBI_sequence_type" not in df.columns:
    raise ValueError(
        "A coluna 'NCBI_sequence_type' não foi encontrada na aba 'tudo'."
    )

# Normalizar os valores da coluna
tipos = (
    df["NCBI_sequence_type"]
    .astype("string")
    .str.strip()
    .str.lower()
)

# Separar as cepas
cromossomos = df.loc[tipos == "chromosome"].copy()
plasmidios = df.loc[tipos == "plasmid"].copy()

# Salvar os cromossomos em um arquivo tabular
cromossomos.to_csv(
    arquivo_cromossomos,
    sep="\t",
    index=False,
    encoding="utf-8-sig"
)

# Salvar os plasmídios em outro arquivo tabular
plasmidios.to_csv(
    arquivo_plasmidios,
    sep="\t",
    index=False,
    encoding="utf-8-sig"
)

# Exibir os resultados
print(f"Arquivo de cepas cromossomais criado: {arquivo_cromossomos}")
print(f"Arquivo de cepas plasmidiais criado: {arquivo_plasmidios}")
print(f"Cepas cromossomais: {len(cromossomos)}")
print(f"Cepas plasmidiais: {len(plasmidios)}")

