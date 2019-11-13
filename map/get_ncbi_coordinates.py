import re

import pandas as pd
from urllib import request
import xml.etree.ElementTree as xml


def get_assembles():
    data = pd.read_csv(r"C:\Users\artem\PycharmProjects\metavirome\assembly_stats.txt")
    data = data[data["Organism"].str.contains("marine m") |
                data["Organism"].str.contains("sea") |
                data["Organism"].str.contains("soil")]
    assembly_url = r"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=assembly&RetMax=25000&id="
    biosample_url = r"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=biosample&RetMax=25000&id="

    with open(f"coordinates.csv", "w") as csv:
        csv.write("id,id2,latitude,longitude\n")

        for i in range(len(data) // 100 + 1):
            current = data[i * 100: (i + 1) * 100]
            ids = _get_biosample_ids(assembly_url, current["id"].astype(str))
            result = request.urlopen(biosample_url + ",".join(ids)).read().decode()
            root = xml.fromstring(result)
            try:
                for DocumentSummary in root.iter("DocumentSummary"):
                    DocumentSummary = _get_biosample_data(DocumentSummary)
                    id = re.search(r' id="(\d+)"', DocumentSummary).group(1)
                    id2 = re.search(r'<Id db="BioSample"(?: is_primary="1")?>([A-Z0-9]+)</Id>', DocumentSummary).group(1)
                    latitude = _get_latitude(DocumentSummary)
                    longitude = _get_longitude(DocumentSummary)
                    csv.write(f"{id},{id2},{latitude},{longitude}\n")
            except Exception as e:
                print(e, id, DocumentSummary)


def _get_biosample_ids(assembly_url, ids):
    result = request.urlopen(assembly_url + ",".join(ids)).read().decode()
    root = xml.fromstring(result)
    ids = [c.text for c in root.iter("BioSampleId")]
    return ids


def _get_biosample_data(data):
    data = [d for d in data.iter("SampleData")]
    assert len(data) == 1, data
    data = data[0].text
    return data


def _get_latitude(data):
    latitude = re.search(r'<Attribute[^/]+latitude[^/]+">-?(\d+(\.\d+)?(.[NS])?)</Attribute>', data)
    if not latitude:
        latitude = re.search(r'(\d+\.\d+ [NS]) \d+\.\d+ [EW]', data)
    if latitude:
        latitude = latitude.group(1)
        if latitude.endswith("N"):
            latitude = latitude[:len(latitude) - 2]
        elif latitude.endswith("S"):
            latitude = "-" + latitude[:len(latitude) - 2]
    return latitude


def _get_longitude(data):
    longitude = re.search(r'<Attribute[^/]+longitude[^/]+">-?(\d+(\.\d+)?(.[WE])?)</Attribute>', data)
    if not longitude:
        longitude = re.search(r'\d+\.\d+ [NS] (\d+\.\d+ [EW])', data)
    if longitude:
        longitude = longitude.group(1)
        if longitude.endswith("E"):
            longitude = longitude[:len(longitude) - 2]
        elif longitude.endswith("W"):
            longitude = "-" + longitude[:len(longitude) - 2]
    return longitude


if __name__ == '__main__':
    get_assembles()
