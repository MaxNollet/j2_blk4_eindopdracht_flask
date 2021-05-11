from dataclasses import dataclass
from gaps.models import Gene
import re


def reader(file, headers):
    """
    Reads the GenePanel file, gets the headers from the first row
    :param file: GenPanelOverzicht | tsv file | first row needs to contain
    the names of the columns
    :param headers: names of the columns to get the data from
    :return:
    """
    data = {}  # is this the correct data structure?
    # Gene(id="")
    # Gene object contains: GeneID
    data_test = []
    with open(file, mode="r", encoding="utf-8-sig") as f:
        # https://stackoverflow.com/questions/17912307/u-ufeff-in-python-string
        # data = [i for i, j in enumerate(f.readline().strip().split(";")) if "GenePanels" in j]
        # how to replace "GenePanels" with h? without setting the entire line in the loop
        # for i, j in enumerate(f.readline().strip().split(";")):
        for i, j in enumerate(f.readline().strip().split("\t")):
            # i is the index and j is the value

            for h in headers:
                # print(h, j)
                if h in j:
                    # data.update({i: j})  # hier staat er \ufeff voor

                    data[i] = [j]
        for line in f:
            gene = Gene(in_genepanel=True)
            for key_index, value in data.items():  # beter naam voor
                # value nodig beetje karig het is nu gwn een list waar alle data van een kolom in komt
                # print(key, value)
                # value.append(line.strip().split("\t")[key_index])
                # if re.findall("(?<=_).+?(?=\])", headers[key_index]) == "ncbi":
                if re.search("NCBI", headers[key_index]):
                    gene.ncbi_gene_id = line.strip().split("\t")[key_index]
                if re.search("HGNC", headers[key_index]):
                    gene.hgnc_symbol = line.strip().split("\t")[key_index]
                value.append(gene)
            data_test.append(gene)
            # print(line)

    # if headers[key_index] == "GeneID_NCBI":
    # value.append(Gline.strip().split("\t")[key_index])
    # value.append(line.strip().split("\t")[keyword])
    # print(line.strip().split("\t")[key_index])
    # add column to selected keyword

    print(len(data.get(1)))
    print(len(data_test))
    print(data_test)

    # print(data.get(0)[1])
    # t = [Gene(ncbi_gene_id=data.get(0)[1])]
    # a = [Gene(ncbi_gene_id=data.get(0)[1])]
    # a[0].hgnc_symbol=data.get(1)[1]
    # print(a)
    # print(a[0])

    # tt = ParseItems(ncbi_gene_id=int(data.get(0)[1]))
    # print(tt)

    # print(len(data.get(0)))
    # print(data)
    for keyword, value in data.items():
        print(keyword)

    print("Voltooid")

@dataclass()
class ParseItems:
    ncbi_gene_id: int
    hgnc_symbol = str
    alias = list
    in_panel = bool


def main():
    # headers = ["GeneID_NCBI", "Symbol_HGNC", "Aliases"]
    headers = ["GeneID_NCBI", "Symbol_HGNC"]
    file = "/Users/lean/Documenten/School/Flask/Course8_project/Parser/GenPanelOverzicht_edited.csv"
    file = "/Users/lean/Documenten/School/Flask/Course8_project/Parser/GenPanelOverzicht_DG-3.1.0_HAN_original_tsv.txt"
    reader(file, headers)
    # f = open(file)
    # print(f.read())


main()
