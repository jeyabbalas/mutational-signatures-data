from collections import OrderedDict
import gzip
from pathlib import Path
import shutil
from typing import List

from natsort import natsorted
import numpy as np
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
        data_id = data_id.groupby("chromosome", sort=False) \
            .apply(pd.DataFrame.sort_values, "chromosome_start") \
            .reset_index(drop=True)
        data_id.to_csv(dir_output / f"{donor_id}", sep="\t", index=False)


def convert_ssms_to_mafs(dir_datasets: Path, dir_output: Path) -> None:
    """
    Converts each SSM dataset in a directory into MAF files.

    :param dir_datasets: directory containing SSM datasets.
    :param dir_output: directory to store MAF files.
    """
    filepaths = list(dir_datasets.glob("*.tsv.gz"))

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
        for data in tqdm(response.iter_content(10 * 1024 ** 2)):
            f.write(data)


def gunzip(gzipped_filepath: Path, gunzipped_filepath: Path) -> None:
    """
    Uncompress a gzipped file.

    :param gzipped_filepath: gzip compressed filepath.
    :param gunzipped_filepath: filepath for unzipped file.
    """
    with gzip.open(gzipped_filepath, "rb") as f_src:
        with open(gunzipped_filepath, "wb") as f_dest:
            shutil.copyfileobj(f_src, f_dest, length=10 * 1024 ** 2)


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
        for substitution in substitution_types:
            for base_3 in nucleotide_bases:
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
        sbs_mutational_spectra[context] = [0] * n_records

    return sbs_mutational_spectra


def index_reference_genome(ref_fasta_filepath: Path) -> pyfaidx.Fasta:
    """
    Returns an indexed FASTA file to quickly lookup subsequences in a genome.

    :param ref_fasta_filepath: filepath of the FASTA file of the reference genome.
    :return: an indexed FASTA file
    """
    return pyfaidx.Fasta(ref_fasta_filepath.as_posix())


def read_sbs_maf_file(filepath: str) -> pd.DataFrame:
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
    return ref_fasta[f"chr{row['chromosome']}"][(pointer - 2):(pointer + 1)].seq.upper()


def standardize_trinucleotide(trinucleotide_ref: str) -> str:
    """
    COSMIC signatures define mutations from a pyrimidine allele (C, T) to any
    other base (C>A, C>G, C>T, T>A, T>C, T>G). If a mutation in the MAF file
    is defined from a purine allele (A, G), then we infer the trinucleotide
    context in the complementary sequence, which would be from a pyrimidine
    allele due to purines and pyrimidines complementing each other in a
    double-stranded DNA.

    :param trinucleotide_ref: trinucleotide sequence seen in the reference genome.
    :return: a pyrimidine-centric trinucleotide sequence.
    """
    complement_seq = {
        'A': 'T',
        'C': 'G',
        'T': 'A',
        'G': 'C'
    }
    purines = ["A", "G"]
    if trinucleotide_ref[1] in purines:
        return f"{complement_seq[trinucleotide_ref[2]]}" \
               f"{complement_seq[trinucleotide_ref[1]]}" \
               f"{complement_seq[trinucleotide_ref[0]]}"
    else:
        return trinucleotide_ref


def standardize_substitution(ref_allele: str,
                             mut_allele: str) -> str:
    """
    COSMIC signatures define mutations from a pyrimidine allele (C, T) to any
    other base (C>A, C>G, C>T, T>A, T>C, T>G). If a mutation in the MAF file
    is defined from a reference purine allele (A, G), then we infer the substituted
    base in the complementary sequence, which would be from a pyrimidine
    allele due to purines and pyrimidines complementing each other in a
    double-stranded DNA.

    :param ref_allele: base in the reference genome.
    :param mut_allele: base in the mutated genome
    :return: substitution string from pyrimidine to any other base.
    """
    complement_seq = {
        'A': 'T',
        'C': 'G',
        'T': 'A',
        'G': 'C'
    }
    purines = ["A", "G"]
    if ref_allele in purines:
        return f"{complement_seq[ref_allele]}>{complement_seq[mut_allele]}"
    else:
        return f"{ref_allele}>{mut_allele}"


def add_instance_to_mutational_spectra(maf_df: pd.DataFrame,
                                       mutational_spectra: OrderedDict[str, List[int]],
                                       ref_fasta: pyfaidx.Fasta,
                                       index: int) -> None:
    """
    Parses each row in a MAF dataframe generated from an ICGC SSM dataset and tabulates a
    mutational spectra count matrix in the form of an ordered dictionary.

    :param maf_df: MAF dataframe generated from an ICGC SSM dataset.
    :param mutational_spectra: an ordered dictionary to tabulat the mutational spectra matrix.
    :param ref_fasta: an indexed reference genome.
    :param index: row index in the mutational spectra matrix to tabulate in the counts.
    """
    nucleotide_bases = ["A", "C", "G", "T"]
    pyrimidine = ["C", "T"]

    for _, row in maf_df.iterrows():
        if ((row["chromosome_start"] != row["chromosome_end"]) or
                (row["reference_genome_allele"] not in nucleotide_bases) or
                (row["mutated_to_allele"] not in nucleotide_bases)):
            continue
        trinucleotide_ref = standardize_trinucleotide(
            get_trinucleotide_ref_from_fasta(row, ref_fasta))
        substitution = standardize_substitution(row["reference_genome_allele"],
                                                row["mutated_to_allele"])

        # sanity checks
        try:
            assert (trinucleotide_ref is not None)
            assert (trinucleotide_ref[1] == substitution[0])
            assert (trinucleotide_ref[1] in pyrimidine)
            assert (substitution[0] in pyrimidine)
        except AssertionError:
            print(f"MAF row: {row['chromosome']}, "
                  f"{row['chromosome_start']}, "
                  f"{row['chromosome_end']}, "
                  f"{row['reference_genome_allele']}, "
                  f"{row['mutated_to_allele']}")
            print(f"FASTA context: {get_trinucleotide_ref_from_fasta(row, ref_fasta)}")
            print(f"Pyrimidine-centric context: {trinucleotide_ref}")
            raise

        mutational_spectra[f"{trinucleotide_ref[0]}[{substitution}]{trinucleotide_ref[2]}"][index] += 1


def write_mutational_spectra(mutational_spectra: OrderedDict,
                             sample_names: List[str],
                             filepath: Path) -> None:
    """
    Writes the mutational spectra matrix data, stored in an ordered dictionary, to a CSV file.

    :param mutational_spectra: mutational spectra matrix data stored in an ordered dictionary.
    :param sample_names: a list of names of the samples.
    :param filepath: name of the CSV file to save the data.
    """
    data = np.stack([np.array(mutational_spectra[substitution]) for substitution in mutational_spectra.keys()])
    index = pd.Series(
        data=mutational_spectra.keys(),
        name="Mutation Types"
    )
    mutational_spectra_df = pd.DataFrame(
        data=data,
        index=index,
        columns=sample_names,
        dtype=int,
    )
    mutational_spectra_df.to_csv(filepath, sep=",", index=True)


def convert_mafs_to_sbs_mutational_spectra(dir_mafs: Path,
                                           ref_fasta_filepath: Path,
                                           filepath_output: Path) -> None:
    """
    Converts all MAF files (one file per sample) in a directory into a mutational spectra
    matrix and saves it as a CSV file.

    :param dir_mafs: a directory containing MAF files.
    :param ref_fasta_filepath: filepath to the reference genome FASTA file.
    :param filepath_output: file path to save the mutational spectra CSV file.
    """
    maf_filepaths = natsorted(list(dir_mafs.glob("*")))
    n_samples = len(maf_filepaths)
    mutational_spectra = init_sbs_mutational_spectra(n_samples)
    ref_fasta = index_reference_genome(ref_fasta_filepath)
    donors = list()

    donor_index = 0
    for maf_filepath in maf_filepaths:
        data_maf = read_sbs_maf_file(maf_filepath)
        add_instance_to_mutational_spectra(data_maf, mutational_spectra, ref_fasta, donor_index)
        donors.append(maf_filepath.name)
        donor_index += 1
    write_mutational_spectra(mutational_spectra, donors, filepath_output)


def convert_maf_dirs_to_sbs_mutational_spectra(dir_maf_dirs: Path,
                                               ref_fasta_filepath: Path,
                                               dir_output: Path) -> None:
    """
    For each directory within the specified directory, this method iterates through all
    MAF files and creates a mutational spectra matrix and saves them as a CSV file.

    :param dir_maf_dirs: a directory of directories, each containing a set of MAF files.
    :param ref_fasta_filepath: filepath to the reference genome FASTA file.
    :param dir_output: directory to save the mutational spectra CSV files.
    """
    for dir_mafs in dir_maf_dirs.iterdir():
        if dir_mafs.is_dir():
            filepath_output = dir_output / f"{dir_mafs.name}.csv"
            convert_mafs_to_sbs_mutational_spectra(dir_mafs, ref_fasta_filepath, filepath_output)
