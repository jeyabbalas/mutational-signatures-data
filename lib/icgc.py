from pathlib import Path

import requests
from tqdm import tqdm


def get_filesize(pql_query: str,
                 datatype: str = "ssm") -> int:
    """
    Calls an ICGC Data Portal API to retrieve the file size of the dataset
    specified by a PQL query and a data type.

    :param pql_query: PQL query to retrieve the dataset of interest.
    :param datatype: data types e.g., "ssm" for simple somatic mutation,
        "donor" for clinical dataset, "cnsm" for copy number somatic mutation,
        etc.
    :return: size of the specified dataset in bytes.
    """
    url = f"https://dcc.icgc.org/api/v1/download/sizePQL?pql={pql_query}"

    response = requests.get(url)
    if response.status_code != 200:
        raise IOError(f"GET {url} resulted in status code {response.status_code}")

    file_sizes = response.json()["fileSize"]
    for dataset in file_sizes:
        if dataset["label"] == datatype:
            return dataset["sizes"]

    raise ValueError(f"GET {url} does not contain the {datatype} data type.")


def get_download_id(pql_query: str,
                    datatype: str = "ssm",
                    output_format: str = "TSV") -> str:
    """
    Calls an ICGC Data Portal API to retrieve a download ID for the dataset
    specified by a PQL query, a data type, and an output format.

    :param pql_query: PQL query to retrieve the dataset of interest.
    :param datatype: data types e.g., "ssm" for simple somatic mutation,
        "donor" for clinical dataset, "cnsm" for copy number somatic mutation,
        etc.
    :param output_format: output data format. Supported formats: ["json", "TSV"].
    :return: a download ID
    """
    info = f"[{{\"key\":\"{datatype}\", \"value\":\"{output_format}\"}}]"
    url = f"https://dcc.icgc.org/api/v1/download/submitPQL?pql={pql_query}&info={info}"

    response = requests.get(url)
    if response.status_code != 200:
        raise IOError(f"GET {url} resulted in status code {response.status_code}")

    return response.json()["downloadId"]


def download_data(output_filepath: Path,
                  download_id: str,
                  file_size: int) -> None:
    """
    Calls an ICGC Data Portal API to download a gzipped file for the dataset
    specified by a download ID.

    :param output_filepath: output file directory
    :param download_id: download ID obtained from API call from get_download_id()
    :param file_size: dataset file size in bytes
    """
    url = f"https://dcc.icgc.org/api/v1/download/{download_id}"
    headers = {"Accept": "application/x-gzip"}
    progress_bar = tqdm(total=file_size, unit="iB", unit_scale=True)

    response = requests.get(url, headers=headers,
                            verify=False, stream=True)
    if response.status_code != 200:
        raise IOError(f"GET {url} resulted in status code {response.status_code}")

    with open(output_filepath.with_suffix(".tsv.gz"), "wb") as f:
        for data in response.iter_content(1024 ** 2):
            progress_bar.update(len(data))
            f.write(data)
    progress_bar.close()


def download_icgc_datasets(output_dir: Path,
                           projects: list[str],
                           datatype: str = "ssm",
                           analysis_type: str = "WGS",
                           output_format: str = "TSV") -> None:
    """
    Download BRCA project datasets from ICGC Data Portal.

    :param output_dir: output directory to download data in.
    :param projects: a list of projects in ICGC to extract data from.
    :param datatype: data types e.g., "ssm" for simple somatic mutation,
        "donor" for clinical dataset, "cnsm" for copy number somatic mutation,
        etc.
    :param analysis_type: data analysis type. E.g., WGS for whole genome sequencing,
        WXS for whole exome sequencing, etc.
    :param output_format: output data format. Supported formats: ["json", "TSV"].
    """
    supported_formats = ["TSV", "json"]
    if output_format not in supported_formats:
        raise ValueError(f"Output format {output_format} isn't supported. "
                         f"Supported formats: {supported_formats}")

    if not output_dir.exists():
        output_dir.mkdir()

    for project in projects:
        pql_query = f"select(*),in(donor.projectId,'{project}')," \
                    f"in(donor.availableDataTypes,'{datatype}')," \
                    f"in(donor.analysisTypes,'{analysis_type}')"

        file_size = get_filesize(pql_query, datatype)
        print(f"Downloading {datatype} data ({(file_size / 1024 ** 2):.2f} MBs) "
              f"from project {project}.")
        download_id = get_download_id(pql_query, datatype, output_format)
        output_filepath = output_dir / f"{project}_{datatype}_{analysis_type}"
        download_data(output_filepath, download_id, file_size)

    print("Done.")
