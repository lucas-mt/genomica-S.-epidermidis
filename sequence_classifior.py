#!/usr/bin/env python3
"""
Consulta os accessions da coluna 'accession' da aba 'tudo' de uma planilha ODS
no NCBI e grava somente os resultados em um arquivo TSV.

Uso:
    python classificar_genbank_ncbi.py "planilha.ods"

Dependências:
    pip install pandas odfpy requests biopython
"""
import sys
import time
import re
from pathlib import Path
from io import StringIO

import pandas as pd
import requests
from Bio import SeqIO

NCBI_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
HEADERS = {"User-Agent": "GenBankAccessionClassifier/1.1 (academic research)"}


def classify_record(record):
    """Classifica usando o registro GenBank retornado pelo NCBI."""
    definition = (record.description or "").strip()
    definition_l = definition.lower()
    keywords = " ".join(str(x) for x in record.annotations.get("keywords", [])).lower()

    plasmid_values = []
    chromosome_values = []
    for feature in record.features:
        if feature.type == "source":
            plasmid_values.extend(feature.qualifiers.get("plasmid", []))
            chromosome_values.extend(feature.qualifiers.get("chromosome", []))

    # Evidências explícitas do registro têm prioridade.
    if plasmid_values or re.search(r"\bplasmid\b", definition_l) or "plasmid" in keywords:
        return "Plasmid"
    if chromosome_values or re.search(r"\bchromosome\b", definition_l):
        return "Chromosome"

    # Para S. epidermidis, registros cromossômicos completos geralmente são > 1 Mb.
    # Se não há evidência de plasmídeo e o registro é grande, classifica como cromossomo provável.
    length = len(record.seq)
    if length >= 1_000_000:
        return "Chromosome (inferred by sequence length)"

    # Não força uma classificação quando os dados do NCBI não bastam.
    if "complete genome" in definition_l or "genome assembly" in definition_l:
        return "Genome/assembly - review"
    if any(term in definition_l for term in ("contig", "scaffold", "whole genome shotgun", "wgs")):
        return "Draft genome/contig - review"
    return "Uncertain - manual review"


def fetch_accession(accession, session):
    params = {
        "db": "nuccore",
        "id": accession,
        "rettype": "gb",
        "retmode": "text",
        "tool": "GenBankAccessionClassifier",
    }
    for attempt in range(4):
        try:
            response = session.get(NCBI_EFETCH, params=params, headers=HEADERS, timeout=45)
            response.raise_for_status()
            text = response.text
            if not text.strip() or text.lstrip().startswith("<"):
                raise ValueError("NCBI não retornou um registro GenBank em texto")
            records = list(SeqIO.parse(StringIO(text), "genbank"))
            if not records:
                raise ValueError("Registro não encontrado ou formato inesperado")
            record = records[0]
            return {
                "NCBI_sequence_type": classify_record(record),
                "NCBI_definition": record.description,
                "NCBI_record_url": f"https://www.ncbi.nlm.nih.gov/nuccore/{accession}",
                "NCBI_check_status": "OK",
            }
        except Exception as exc:
            if attempt == 3:
                return {
                    "NCBI_sequence_type": "Not classified",
                    "NCBI_definition": "",
                    "NCBI_record_url": f"https://www.ncbi.nlm.nih.gov/nuccore/{accession}",
                    "NCBI_check_status": f"ERROR: {str(exc)[:180]}",
                }
            time.sleep(2 ** attempt)
    return {}


def main():
    if len(sys.argv) < 2:
        print('Uso: python classificar_genbank_ncbi.py "arquivo.ods"')
        raise SystemExit(2)

    input_path = Path(sys.argv[1]).expanduser().resolve()
    if not input_path.exists():
        raise FileNotFoundError(f"Arquivo não encontrado: {input_path}")

    # Lê SOMENTE a aba 'tudo', conforme solicitado.
    try:
        df = pd.read_excel(input_path, sheet_name="tudo", engine="odf")
    except ValueError as exc:
        raise ValueError("A planilha não possui uma aba chamada 'tudo'. Confira o nome da aba.") from exc

    # Localiza a coluna accession sem diferenciar maiúsculas/minúsculas ou espaços.
    accession_col = next(
        (col for col in df.columns if str(col).strip().lower() in
         ("accession", "genbank accession", "genbank_accession")),
        None,
    )
    if accession_col is None:
        raise ValueError(
            "Não encontrei a coluna 'accession' na aba 'tudo'. "
            f"Colunas encontradas: {', '.join(map(str, df.columns))}"
        )

    accessions = []
    seen = set()
    for raw in df[accession_col]:
        accession = "" if pd.isna(raw) else str(raw).strip()
        if accession and accession.lower() not in {"nan", "none"} and accession not in seen:
            accessions.append(accession)
            seen.add(accession)

    if not accessions:
        raise ValueError("Nenhum accession preenchido foi encontrado na coluna da aba 'tudo'.")

    session = requests.Session()
    results = []
    total = len(accessions)
    for i, accession in enumerate(accessions, start=1):
        result = fetch_accession(accession, session)
        results.append({"accession": accession, **result})
        time.sleep(0.36)  # limite de frequência do NCBI sem API key
        if i % 25 == 0 or i == total:
            print(f"Consultados {i}/{total} accessions únicos da aba 'tudo'")

    result_df = pd.DataFrame(results, columns=[
        "accession", "NCBI_sequence_type", "NCBI_definition",
        "NCBI_record_url", "NCBI_check_status",
    ])
    output_path = input_path.with_name(input_path.stem + "_informacoes_NCBI.tsv")
    result_df.to_csv(output_path, sep="\t", index=False, encoding="utf-8")

    print(f"\nConcluído. Resultados salvos em: {output_path}")
    print("A saída contém apenas os accessions únicos encontrados na aba 'tudo' e os resultados da consulta ao NCBI.")
    print("Confira manualmente os registros marcados como 'review' ou 'Uncertain'.")


if __name__ == "__main__":
    main()
