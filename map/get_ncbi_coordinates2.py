from urllib import request
import xml.etree.ElementTree as xml
import re


def _get_ids(file_name, id_tag):
    root = xml.parse(file_name).getroot()
    ids = ",".join(c.text for c in root.iter(id_tag))
    return ids


def get_assemble_ids():
    ids = _get_ids("esearch.xml", "Id")
    return ids


def get_biosample_ids():
    ids = _get_ids("assembly.xml", "BioSampleId")
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
    term = "sea"
    file_name = "sea"
    esearch_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=assembly&RetMax=25000&term={term}"
    assembly_url = r"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=assembly&RetMax=25000&id="
    biosample_url = r"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=biosample&RetMax=25000&id="

    request.urlretrieve(esearch_url, "esearch.xml")
    ids = get_assemble_ids()

    request.urlretrieve(assembly_url + ids, f"assembly.xml")
    ids = get_biosample_ids()

    request.urlretrieve(biosample_url + ids, "biosample_ansi.xml")
    with open("biosample.xml", "w", encoding="utf-8") as reencoded:
        with open("biosample_ansi.xml", "r", encoding="ansi") as ncbi_metadata_is_a_mess:
            reencoded.writelines(ncbi_metadata_is_a_mess.readlines())

    root = xml.parse("biosample.xml").getroot()

    with open(f"coordinates_{file_name}.csv", "w") as csv:
        try:
            csv.write("id,id2,latitude,longitude\n")
            for data in root.iter("DocumentSummary"):
                data = _get_biosample_data(data)
                id = re.search(r' id="(\d+)"', data).group(1)
                id2 = re.search(r'<Id db="BioSample"(?: is_primary="1")?>([A-Z0-9]+)</Id>', data).group(1)
                latitude = _get_latitude(data)
                longitude = _get_longitude(data)
                csv.write(f"{id},{id2},{latitude},{longitude}\n")
        except Exception as e:
            print(e, id, data)

