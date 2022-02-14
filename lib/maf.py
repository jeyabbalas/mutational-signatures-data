from collections import OrderedDict
import gzip
from pathlib import Path
import shutil
from typing import List

import pandas as pd
import pyfaidx
import requests
from tqdm import tqdm


def read_ssm_dataset(filepath: Path) -> pd.DataFrame:
    """
    Reads an ICGC SSM file as a pandas dataframe.

    :param filepath: file path to the SSM dataset.
    :return: pandas dataframe selecting only the columns relevant to mutational
        signatures analysis.
    """
    select_columns = [
        "icgc_mutation_id", "project_code", "icgc_donor_id",
        "chromosome", "chromosome_start", "chromosome_end",
        "assembly_version", "mutation_type", "reference_genome_allele",
        "mutated_to_allele",
    ]

    return pd.read_csv(filepath, usecols=select_columns, sep="\t")


def clean_ssm_dataset(data: pd.DataFrame) -> pd.DataFrame:
    """
    Keeps only one variant per donor ID and drops the rest. The repeats are due to
    the SnpEff annotation tool, which is initially irrelevant for signature analysis.

    :param data: a dataframe of SSM.
    :return: a dataframe of SSM without repeats
    """
    return data.drop_duplicates(subset=["icgc_donor_id", "icgc_mutation_id"]).reset_index(drop=True)


def segregate_ids_and_save_as_maf(data: pd.DataFrame,
                                  dir_output: Path) -> None:
    """
    Takes an ICGC SSM dataset, groups them by donor ID, then for each donor ID,
    sorts the records by chromosome number and then by chromosome start position,
    and finally writes this dataset as an MAF file.

    :param data: SSM dataframe
    :param dir_output: output directory for the MAF files.
    """
    for donor_id in data["icgc_donor_id"].unique():
        data_id = data.loc[data["icgc_donor_id"] == donor_id]
        data_id = data_id.loc[pd.to_numeric(data_id["chromosome"], errors="coerce").sort_values().index]
        data_id = data_id.groupby("chromosome", sort=False)\
            .apply(pd.DataFrame.sort_values, "chromosome_start")\
            .reset_index(drop=True)
        data_id.to_csv(dir_output / f"{donor_id}", sep="\t", index=False)


def convert_ssms_to_mafs(dir_datasets: Path, dir_output: Path) -> None:
    """
    Converts each SSM dataset in a directory into MAF files.

    :param dir_datasets: directory containing SSM datasets.
    :param dir_output: directory to store MAF files.
    """
    filepaths = list(dir_datasets.glob("*.tsv.gz"))

    dir_output = dir_output / (dir_datasets.name + "_MAF")
    if not dir_output.exists():
        dir_output.mkdir()

    for filepath in filepaths:
        data = read_ssm_dataset(filepath)
        data = clean_ssm_dataset(data)
        dir_output_file = dir_output / (filepath.stem.split("_")[0])
        if not dir_output_file.exists():
            dir_output_file.mkdir()
        segregate_ids_and_save_as_maf(data, dir_output_file)


def download_grch37(filepath: Path) -> None:
    """
    Downloads a compressed FASTA file of the reference genome GRCh37 from the
    UCSC Genome Browser API.

    :param filepath: output directory.
    """
    url = "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz"
    headers = {"Accept": "application/x-gzip"}

    response = requests.get(url, headers=headers,
                            verify=False, stream=True)
    if response.status_code != 200:
        raise IOError(f"GET {url} resulted in status code {response.status_code}")

    with open(filepath, "wb") as f:
        for data in tqdm(response.iter_content(10*1024**2)):
            f.write(data)


def gunzip(gzipped_filepath: Path, gunzipped_filepath: Path) -> None:
    """
    Uncompress a gzipped file.

    :param gzipped_filepath: gzip compressed filepath.
    :param gunzipped_filepath: filepath for unzipped file.
    """
    with gzip.open(gzipped_filepath, "rb") as f_src:
        with open(gunzipped_filepath, "wb") as f_dest:
            shutil.copyfileobj(f_src, f_dest, length=10*1024**2)


def get_sbs_trinucleotide_contexts() -> List[str]:
    """
    Returns a list of trinucleotide context for single base substitutions (SBS)
    for constructing a COSMIC mutational spectra matrix.

    :return: a list of SBS trinucleotide contexts.
    """
    sbs_trinucleotide_contexts = []
    nucleotide_bases = ["A", "C", "G", "T"]
    substitution_types = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]

    for base_5 in nucleotide_bases:
        for base_3 in nucleotide_bases:
            for substitution in substitution_types:
                sbs_trinucleotide_contexts.append(f"{base_5}[{substitution}]{base_3}")

    return sbs_trinucleotide_contexts


def init_sbs_mutational_spectra(n_records: int) -> OrderedDict[str, List[int]]:
    """
    Initilizes an ordered dictionary with SBS trinucleotide context as keys and
    a list of counts, one for each sample.

    :param n_records: number of samples to record in the mutational spectra matrix.
    :return: an ordered dictionary of trinucleotide context and a list of counts
        initialized to zeros.
    """
    sbs_mutational_spectra = OrderedDict()
    sbs_trinucleotide_contexts = get_sbs_trinucleotide_contexts()

    for context in sbs_trinucleotide_contexts:
        sbs_mutational_spectra[context] = [0]*n_records

    return sbs_mutational_spectra


def index_reference_genome(ref_fasta_filepath: Path) -> pyfaidx.Fasta:
    """
    Returns an indexed FASTA file to quickly lookup subsequences in a genome.

    :param ref_fasta_filepath: filepath of the FASTA file of the reference genome.
    :return: an indexed FASTA file
    """
    return pyfaidx.Fasta(ref_fasta_filepath)


def read_sbs_maf_file(filepath: Path) -> pd.DataFrame:
    """
    Reads only single base substitutions from an MAF file generated from an
    ICGC SSM dataset.

    :param filepath: file path to the MAF file.
    :return: Pandas dataframe with only single base substitutions.
    """
    data = pd.read_csv(filepath, sep="\t")
    data = data.loc[data["mutation_type"] == "single base substitution"].reset_index(drop=True)
    return data


def get_trinucleotide_ref_from_fasta(row: pd.Series,
                                     ref_fasta: pyfaidx.Fasta) -> str:
    """
    Returns the trinucleotides (5' base, reference allele, 3' base) around the
    mutation described by the row.

    :param row: a pandas row of the MAF file.
    :param ref_fasta: an indexed FASTA of the reference genome.
    :return: trinucleotide context for the mutation described by the row.
    """
    pointer = int(row["chromosome_start"])
    """
    '-2' and not '-1' because genomes are indexed starting from 1 but Python data
    structures are indexed starting from 0.
    """
    return ref_fasta[f"chr_{row['chromosome']}"][(pointer-2):(pointer+1)].seq.upper()