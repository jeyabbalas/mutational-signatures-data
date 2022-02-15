from pathlib import Path

import pandas as pd

from maf import read_sbs_maf_file


def is_mutation_in_gene_panel(panel_df,
                              chromosome: str,
                              position_start: int,
                              position_end: int) -> bool:
    """
    Checks if the defined mutation is measured by the gene panel.

    :param panel_df: gene panel BED file.
    :param chromosome: chromosome name of the mutation.
    :param position_start: starting position of the mutation.
    :param position_end: end position of the mutation.
    :return: True if the mutation region is measured in the gene panel.
        False otherwise.
    """
    panel_matches = panel_df.loc[((panel_df["Chromosome"] == chromosome) &
                                  (panel_df["Start_Position"] <= position_start) &
                                  (panel_df["End_Position"] >= position_end))]
    return len(panel_matches) > 0


def return_panel_rows(panel_df: pd.DataFrame,
                      wg_maf: pd.DataFrame) -> pd.DataFrame:
    """
    Takes an MAF dataframe from WGS analysis and downsamples it to only return the
    mutations covered by the gene panel (as described by its BED file).

    :param panel_df: gene panel's BED file.
    :param wg_maf: MAF from WGS analysis.
    :return:
    """
    panel_rows = wg_maf.apply(lambda row: is_mutation_in_gene_panel(panel_df,
                                                                    row["chromosome"],
                                                                    row["chromosome_start"],
                                                                    row["chromosome_end"]),
                              axis=1)
    return wg_maf.loc[panel_rows]


def convert_wg_mafs_to_panel_mafs(dir_wg_mafs: Path,
                                  panel_df: pd.DataFrame,
                                  dir_panel_mafs: Path) -> None:
    """
    Converts each MAF file, generated from WGS analysis, in a directory and
    downsamples it to only MAF files containing only the mutations covered
    by the gene panel (as described by its BED file).

    :param dir_wg_mafs: a directory containing a set of MAF files from WGS
        analysis.
    :param panel_df: BED file of the gene panel.
    :param dir_panel_mafs: output directory.
    """
    maf_filepaths = list(dir_wg_mafs.glob("*"))

    for maf_filepath in maf_filepaths:
        filepath_panel = dir_panel_mafs / f"{maf_filepath.stem}"
        data_maf = read_sbs_maf_file(maf_filepath)
        data_maf = return_panel_rows(panel_df, data_maf)
        if len(data_maf) > 0:
            data_maf.to_csv(filepath_panel, sep="\t", index=False)


def panel_maf_dirs_from_wgs_maf_dirs(dir_wg_maf_dirs: Path,
                                     panel_df: pd.DataFrame,
                                     dir_panel_maf_dirs: Path) -> None:
    """
    For each directory within the specified directory, this method iterates through
    all MAF files, generated from WGS analysis, and downsamples it to create MAF
    files containing mutations covered by the gene panel (as described by its BED
    file).

    :param dir_wg_maf_dirs: a directory of directories, each containing a set of MAF files
        from WGS analysis.
    :param panel_df: BED file of the gene panel.
    :param dir_panel_maf_dirs: output directory.
    """
    for dir_wg_mafs in dir_wg_maf_dirs.iterdir():
        dir_panel_mafs = dir_panel_maf_dirs / dir_wg_mafs.name
        if not dir_panel_mafs.exists():
            dir_panel_mafs.mkdir()
        convert_wg_mafs_to_panel_mafs(dir_wg_mafs, panel_df, dir_panel_mafs)
