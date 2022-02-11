from pathlib import Path

from lib import icgc, maf


if __name__ == '__main__':
    dir_data = Path.cwd().parent / "data"
    dir_wgs = dir_data / "WGS"

    projects = ["BRCA-EU", "BRCA-FR", "BRCA-UK", "BRCA-US"]
    datatype = "ssm"
    analysis_type = "WGS"
    output_format = "TSV"

    if not dir_data.exists():
        dir_data.mkdir()
    if not dir_wgs.exists():
        dir_wgs.mkdir()

    icgc.download_icgc_datasets(dir_wgs, projects, datatype, analysis_type, output_format)
    maf.convert_ssms_to_mafs(dir_wgs, dir_data)
