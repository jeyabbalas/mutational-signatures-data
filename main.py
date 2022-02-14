from pathlib import Path

from lib import icgc, maf


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

    maf.convert_ssms_to_mafs(dir_wgs, dir_data)
    fpath_compressed_grch37 = dir_data / "hg19.fa.gz"
    maf.download_grch37(fpath_compressed_grch37)
    fpath_grch37 = dir_data / "hg19.fa"
    maf.gunzip(fpath_compressed_grch37, fpath_grch37)
