from pathlib import Path

import pandas as pd

from lib import icgc, maf, msk_impact_410


if __name__ == '__main__':
    dir_data = Path.cwd().parent / "data"
    if not dir_data.exists():
        dir_data.mkdir()

    projects = ["BRCA-EU", "BRCA-FR", "BRCA-UK", "BRCA-US"]
    datatype = "ssm"
    analysis_type = "WGS"
    output_format = "TSV"

    dir_wgs = dir_data / "WGS"
    icgc.download_icgc_datasets(dir_wgs, projects, datatype, analysis_type, output_format)

    dir_maf_dirs = dir_data / "WGS_MAFs"
    if not dir_maf_dirs.exists():
        dir_maf_dirs.mkdir()
    maf.convert_ssms_to_mafs(dir_wgs, dir_maf_dirs)

    fpath_compressed_grch37 = dir_data / "hg19.fa.gz"
    maf.download_grch37(fpath_compressed_grch37)
    fpath_grch37 = dir_data / "hg19.fa"
    maf.gunzip(fpath_compressed_grch37, fpath_grch37)

    dir_spectra = dir_data / "mutational_spectra_wgs"
    if not dir_spectra.exists():
        dir_spectra.mkdir()

    maf.convert_maf_dirs_to_sbs_mutational_spectra(dir_maf_dirs, fpath_grch37, dir_spectra)

    dir_panel_maf_dirs = dir_data / "MSK_IMPACT_410_MAFs"
    if not dir_panel_maf_dirs.exists():
        dir_panel_maf_dirs.mkdir()
    fpath_msk_impact_bed = dir_data / "MSK-IMPACT410.bed"
    msk_impact_df = pd.read_csv(fpath_msk_impact_bed, sep="\t")
    msk_impact_410.panel_maf_dirs_from_wgs_maf_dirs(dir_maf_dirs, msk_impact_df, dir_panel_maf_dirs)

    dir_panel_spectra = dir_data / "mutational_spectra_panel"
    if not dir_panel_spectra.exists():
        dir_panel_spectra.mkdir()
    maf.convert_maf_dirs_to_sbs_mutational_spectra(dir_panel_maf_dirs, fpath_grch37, dir_panel_spectra)
