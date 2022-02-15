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