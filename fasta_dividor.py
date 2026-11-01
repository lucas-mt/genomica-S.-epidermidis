import argparse
import re
import sys
import unicodedata
from pathlib import Path

import pandas as pd


def normalizar(valor):
    """Normaliza um texto para facilitar a comparação."""
    if pd.isna(valor):
        return ""

    texto = str(valor).strip().lower()
    texto = unicodedata.normalize("NFKD", texto)
    texto = "".join(
        c for c in texto
        if not unicodedata.combining(c)
    )

    return re.sub(r"\s+", " ", texto).strip()


def normalizar_accession(valor):
    """Normaliza um accession para comparação."""
    valor = normalizar(valor)

    # Remove espaços e o prefixo accn|, se existir.
    valor = re.sub(r"^>?\s*accn\|", "", valor)
    valor = valor.split("|")[0].strip()

    return valor.upper()


def criar_aliases_accession(accession):
    """
    Cria formas de comparação do accession.
    Também permite comparar accessions com e sem versão,
    por exemplo, CP118842.1 e CP118842.
    """
    accession = normalizar_accession(accession)

    if not accession:
        return set()

    aliases = {accession}

    sem_versao = re.sub(r"\.\d+$", "", accession)
    aliases.add(sem_versao)

    return aliases


def ler_planilha(caminho_planilha):
    """
    Lê a aba 'tudo' e cria índices por accession e cepa.
    A classificação vem exclusivamente de NCBI_sequence_type.
    """
    caminho_planilha = Path(caminho_planilha)

    if not caminho_planilha.exists():
        raise FileNotFoundError(
            f"Planilha não encontrada: {caminho_planilha}"
        )

    df = pd.read_excel(
        caminho_planilha,
        sheet_name="tudo",
        dtype=str
    )

    # Normaliza os nomes das colunas.
    colunas = {
        normalizar(coluna): coluna
        for coluna in df.columns
    }

    obrigatorias = [
        "accession",
        "ncbi_sequence_type",
        "strain"
    ]

    faltantes = [
        coluna
        for coluna in obrigatorias
        if coluna not in colunas
    ]

    if faltantes:
        raise ValueError(
            "As seguintes colunas não foram encontradas "
            f"na aba 'tudo': {', '.join(faltantes)}. "
            f"Colunas disponíveis: {list(df.columns)}"
        )

    col_accession = colunas["accession"]
    col_tipo = colunas["ncbi_sequence_type"]
    col_strain = colunas["strain"]

    # Cada registro representa uma linha da planilha.
    registros = []

    for _, linha in df.iterrows():
        accession = normalizar_accession(
            linha[col_accession]
        )

        tipo_original = (
            str(linha[col_tipo]).strip()
            if pd.notna(linha[col_tipo])
            else ""
        )

        tipo = normalizar(tipo_original)

        strain_original = (
            str(linha[col_strain]).strip()
            if pd.notna(linha[col_strain])
            else ""
        )

        strain = normalizar(strain_original)

        # Classificações com "~" são consideradas errôneas.
        if "~" in tipo_original:
            classificacao = "ignorar"

        elif tipo == "chromosome":
            classificacao = "chromosome"

        elif tipo == "plasmid":
            classificacao = "plasmid"

        else:
            classificacao = "outro"

        registro = {
            "accession": accession,
            "strain_original": strain_original,
            "strain": strain,
            "classificacao": classificacao,
            "tipo_original": tipo_original,
        }

        registros.append(registro)

    # Índice de accessions.
    indice_accession = {}

    for registro in registros:
        accession = registro["accession"]

        if not accession:
            continue

        for alias in criar_aliases_accession(accession):
            indice_accession.setdefault(alias, []).append(
                registro
            )

    # Índice de cepas.
    indice_strain = {}

    for registro in registros:
        strain = registro["strain"]

        if not strain:
            continue

        indice_strain.setdefault(strain, []).append(
            registro
        )

    return indice_accession, indice_strain, len(registros)


def extrair_accession(cabecalho):
    """
    Extrai o accession de cabeçalhos como:

    >accn|CP118842| Staphylococcus epidermidis strain ...
    """
    correspondencia = re.match(
        r"^>accn\|([^|\s]+)",
        cabecalho.strip(),
        flags=re.IGNORECASE
    )

    if correspondencia:
        return normalizar_accession(
            correspondencia.group(1)
        )

    return ""


def corresponde_strain(cabecalho, strain):
    """
    Verifica se a cepa aparece no cabeçalho.

    Usa limites de palavras para evitar, por exemplo,
    que a cepa 704 seja confundida com a cepa 7049.
    """
    strain = normalizar(strain)
    texto = normalizar(cabecalho)

    if not strain:
        return False

    padrao = (
        r"(?<!\w)"
        + re.escape(strain)
        + r"(?!\w)"
    )

    return re.search(padrao, texto) is not None


def identificar_por_strain(cabecalho, indice_strain):
    """
    Tenta identificar um registro pela cepa.

    Se uma cepa estiver associada a mais de um tipo,
    usa o descritor do cabeçalho somente para resolver
    a ambiguidade entre registros daquela mesma cepa.

    A classificação final continua vindo da planilha.
    """
    candidatos = []

    # Procura todas as cepas que aparecem no cabeçalho.
    for strain, registros in indice_strain.items():
        if corresponde_strain(cabecalho, strain):
            candidatos.extend(registros)

    if not candidatos:
        return None, "strain_nao_encontrada"

    # Mantém apenas as cepas mais específicas.
    maior_tamanho = max(
        len(registro["strain"])
        for registro in candidatos
    )

    candidatos = [
        registro
        for registro in candidatos
        if len(registro["strain"]) == maior_tamanho
    ]

    # Remove registros duplicados.
    unicos = {}

    for registro in candidatos:
        chave = (
            registro["accession"],
            registro["strain"],
            registro["classificacao"],
        )
        unicos[chave] = registro

    candidatos = list(unicos.values())

    # Se houver registros errôneos para a cepa,
    # eles não podem ser utilizados para classificação.
    candidatos_validos = [
        registro
        for registro in candidatos
        if registro["classificacao"] != "ignorar"
    ]

    if not candidatos_validos:
        return None, "strain_associada_a_registros_erroneos"

    candidatos = candidatos_validos

    # Se todas as possibilidades possuem a mesma
    # classificação, a identificação é suficiente.
    classificacoes = {
        registro["classificacao"]
        for registro in candidatos
    }

    if len(classificacoes) == 1:
        return candidatos[0], "strain"

    # Caso existam cromossomos e plasmídeos para a
    # mesma cepa, usa o descritor como critério auxiliar.
    # Ele NÃO determina a classificação: apenas seleciona
    # entre registros existentes na planilha.
    cabecalho_normalizado = normalizar(cabecalho)

    if re.search(r"\bchromosome\b", cabecalho_normalizado):
        desejada = "chromosome"

    elif re.search(r"\bplasmid\b", cabecalho_normalizado):
        desejada = "plasmid"

    else:
        return None, "strain_ambiguo"

    candidatos_tipo = [
        registro
        for registro in candidatos
        if registro["classificacao"] == desejada
    ]

    if not candidatos_tipo:
        return None, "strain_sem_correspondencia_de_tipo"

    # Confirma que todos os candidatos restantes
    # possuem a mesma classificação.
    tipos_restantes = {
        registro["classificacao"]
        for registro in candidatos_tipo
    }

    if len(tipos_restantes) != 1:
        return None, "strain_ambiguo"

    return candidatos_tipo[0], "strain"


def classificar_cabecalho(
    cabecalho,
    indice_accession,
    indice_strain
):
    """
    Primeiro tenta localizar o accession.
    Se não conseguir, tenta localizar pela cepa.
    """
    accession = extrair_accession(cabecalho)

    if accession:
        aliases = criar_aliases_accession(accession)

        candidatos = []

        for alias in aliases:
            candidatos.extend(
                indice_accession.get(alias, [])
            )

        # Elimina duplicações.
        unicos = {}

        for registro in candidatos:
            chave = (
                registro["accession"],
                registro["strain"],
                registro["classificacao"],
            )
            unicos[chave] = registro

        candidatos = list(unicos.values())

        if candidatos:
            # Se houver classificação errônea para o
            # accession, ele deve ser ignorado.
            if any(
                registro["classificacao"] == "ignorar"
                for registro in candidatos
            ):
                return "ignorar", "accession_erroneo"

            classificacoes = {
                registro["classificacao"]
                for registro in candidatos
            }

            if len(classificacoes) == 1:
                return (
                    candidatos[0]["classificacao"],
                    "accession"
                )

            # A mesma versão-base pode corresponder
            # a registros diferentes. Tenta desambiguar
            # usando a cepa, sem abandonar a origem
            # da classificação na planilha.
            registro, metodo = identificar_por_strain(
                cabecalho,
                indice_strain
            )

            if registro is not None:
                return registro["classificacao"], metodo

            return "nao_classificado", "accession_ambiguo"

    # Sem accession identificável ou correspondente:
    # tenta localizar a cepa na planilha.
    registro, metodo = identificar_por_strain(
        cabecalho,
        indice_strain
    )

    if registro is not None:
        return registro["classificacao"], metodo

    return "nao_classificado", metodo


def ler_fasta(caminho_fasta):
    """
    Lê o FASTA preservando o cabeçalho e a sequência
    original de cada registro.
    """
    caminho_fasta = Path(caminho_fasta)

    if not caminho_fasta.exists():
        raise FileNotFoundError(
            f"Arquivo FASTA não encontrado: {caminho_fasta}"
        )

    cabecalho = None
    sequencia = []

    with caminho_fasta.open(
        "r",
        encoding="utf-8",
        errors="replace"
    ) as arquivo:

        for numero_linha, linha in enumerate(arquivo, 1):
            linha = linha.rstrip("\r\n")

            if linha.startswith(">"):
                if cabecalho is not None:
                    yield cabecalho, sequencia

                cabecalho = linha
                sequencia = []

            elif cabecalho is not None:
                sequencia.append(linha)

            elif linha.strip():
                raise ValueError(
                    "Foi encontrado texto antes do primeiro "
                    f"cabeçalho FASTA, na linha {numero_linha}."
                )

    if cabecalho is not None:
        yield cabecalho, sequencia


def processar(planilha, fasta):
    caminho_fasta = Path(fasta)

    indice_accession, indice_strain, total_planilha = (
        ler_planilha(planilha)
    )

    # Os arquivos serão criados na mesma pasta do FASTA.
    pasta_saida = caminho_fasta.resolve().parent

    arquivo_cromossomos = pasta_saida / "cromossomos.fasta"
    arquivo_plasmideos = pasta_saida / "plasmideos.fasta"
    arquivo_relatorio = pasta_saida / "relatorio_separacao.tsv"

    contagens = {
        "chromosome": 0,
        "plasmid": 0,
        "ignorar": 0,
        "outro": 0,
        "nao_classificado": 0,
    }

    metodos = {
        "accession": 0,
        "strain": 0,
    }

    # O relatório registra os registros que não foram
    # incluídos nos FASTAs ou que foram ignorados.
    linhas_relatorio = []

    with (
        arquivo_cromossomos.open(
            "w", encoding="utf-8", newline="\n"
        ) as saida_cromossomos,
        arquivo_plasmideos.open(
            "w", encoding="utf-8", newline="\n"
        ) as saida_plasmideos
    ):

        for numero, (cabecalho, sequencia) in enumerate(
            ler_fasta(caminho_fasta), 1
        ):
            classificacao, metodo = classificar_cabecalho(
                cabecalho,
                indice_accession,
                indice_strain
            )

            if metodo in metodos:
                metodos[metodo] += 1

            contagens[classificacao] = (
                contagens.get(classificacao, 0) + 1
            )

            if classificacao == "chromosome":
                destino = saida_cromossomos

            elif classificacao == "plasmid":
                destino = saida_plasmideos

            else:
                destino = None

            if destino is not None:
                destino.write(cabecalho + "\n")

                if sequencia:
                    destino.write("\n".join(sequencia) + "\n")

            else:
                linhas_relatorio.append({
                    "registro_fasta": numero,
                    "accession": extrair_accession(cabecalho),
                    "classificacao": classificacao,
                    "motivo": metodo,
                    "cabecalho": cabecalho,
                })

    # Gera um relatório tabular com registros ignorados,
    # outros tipos de sequência e registros não identificados.
    pd.DataFrame(
        linhas_relatorio,
        columns=[
            "registro_fasta",
            "accession",
            "classificacao",
            "motivo",
            "cabecalho",
        ]
    ).to_csv(
        arquivo_relatorio,
        sep="\t",
        index=False,
        encoding="utf-8"
    )

    print("\nProcessamento concluído!")
    print(f"Registros na planilha: {total_planilha}")
    print(f"Registros classificados por accession: {metodos['accession']}")
    print(f"Registros classificados por cepa: {metodos['strain']}")
    print()
    print(f"Cromossomos: {contagens['chromosome']}")
    print(f"Plasmídeos: {contagens['plasmid']}")
    print(f"Ignorados por erro (~): {contagens['ignorar']}")
    print(f"Outros tipos: {contagens['outro']}")
    print(f"Não classificados: {contagens['nao_classificado']}")
    print()
    print(f"FASTA de cromossomos: {arquivo_cromossomos}")
    print(f"FASTA de plasmídeos: {arquivo_plasmideos}")
    print(f"Relatório: {arquivo_relatorio}")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Separa um arquivo FASTA em cromossomos e "
            "plasmídeos usando a aba 'tudo' de uma planilha."
        )
    )

    parser.add_argument(
        "planilha",
        help="Caminho da planilha .ods ou .xlsx"
    )

    parser.add_argument(
        "fasta",
        help="Caminho do arquivo FASTA original"
    )

    args = parser.parse_args()

    try:
        processar(args.planilha, args.fasta)

    except Exception as erro:
        print(
            f"ERRO: {erro}",
            file=sys.stderr
        )
        sys.exit(1)


if __name__ == "__main__":
    main()

