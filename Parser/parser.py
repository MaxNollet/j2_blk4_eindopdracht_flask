from dataclasses import dataclass
from gaps.models import Gene


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

    with open(file, mode="r", encoding="utf-8-sig") as f:
        # https://stackoverflow.com/questions/17912307/u-ufeff-in-python-string
        # data = [i for i, j in enumerate(f.readline().strip().split(";")) if "GenePanels" in j]
        # how to replace "GenePanels" with h? without setting the entire line in the loop
        # for i, j in enumerate(f.readline().strip().split(";")):
        for i, j in enumerate(f.readline().strip().split("\t")):
            # i is the index and j is the value
            for h in headers:
                if h in j:
                    # data.update({i: j})  # hier staat er \ufeff voor
                    data[i] = [j]
        for line in f:
            for key, value in data.items():
                # print(key, value)
                value.append(line.strip().split("\t")[key])
                # add column to selected keyword
                # print(value)

    print(data)

    print(data.get(0)[1])
    t = [Gene(ncbi_gene_id=data.get(0)[1])]

    # print(len(data.get(0)))
    # print(data)
    for key, value in data.items():
        print(key)

    print("Voltooid")


@dataclass()
class parseItems:
    genid: int
    sym_gene = str
    alias = list
    PanelSymbol = str


def main():
    # headers = ["GeneID_NCBI", "Symbol_HGNC", "Aliases"]
    headers = ["GeneID_NCBI", "Symbol_HGNC"]
    file = "/Users/lean/Documenten/School/Flask/Course8_project/Parser/GenPanelOverzicht_edited.csv"
    file = "/Users/lean/Documenten/School/Flask/Course8_project/Parser/GenPanelOverzicht_DG-3.1.0_HAN_original_tsv.txt"
    reader(file, headers)
    # f = open(file)
    # print(f.read())


main()
