from pathlib import Path

import pandas as pd


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