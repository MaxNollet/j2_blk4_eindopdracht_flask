from dataclasses import dataclass


def reader(file, headers):
    """
    Reads the GenePanel file, gets the headers from the first row
    :param file: GenPanelOverzicht | csv file | first row needs to contain
    the names of the columns
    :param headers: names of the columns to get the data from
    :return:
    """
    data = {}  # test
    data1 = {}  # is this the correct data structure?
    with open(file, mode="r", encoding="utf-8-sig") as f:
        # https://stackoverflow.com/questions/17912307/u-ufeff-in-python-string
        # data = [i for i, j in enumerate(f.readline().strip().split(";")) if "GenePanels" in j]
        # how to replace "GenePanels" with h? without setting the entire line in the loop
        for i, j in enumerate(f.readline().strip().split(";")):
            # i is the index and j is the value
            for h in headers:
                if h in j:
                    data.update({i: j})  # hier staat er \ufeff voor
                    data1[i] = [j]
        for line in f:
            for key, value in data1.items():
                # print(key, value)
                value.append(line.strip().split(";")[key])

    print(data1)
    print(len(data1.get(0)))
    print(data)

    print("Voltooid")


@dataclass()
class parseItems:
    genid: int
    sym_gene = str
    alias = list
    PanelSymbol = str


def main():
    headers = ["GeneID_NCBI", "Symbol_HGNC", "Aliases"]
    file = "/Users/lean/Documenten/School/Flask/Course8_project/Parser/GenPanelOverzicht_edited.csv"
    reader(file, headers)
    # f = open(file)
    # print(f.read())


main()
