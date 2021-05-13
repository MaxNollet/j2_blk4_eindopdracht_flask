from dataclasses import dataclass
from gaps.models import Gene, Alias, FileInfo, Genepanel, GenepanelSymbol
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
    t_data = {}
    t = []  # test
    tt = []  # object test
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
            p_symbol = GenepanelSymbol()
            id = False

            for key_index, value in data.items():  # beter naam voor

                # print(headers[key_index])
                # print(value)
                # value nodig beetje karig het is nu gwn een list waar alle data van een kolom in komt
                # print(key, value)
                # value.append(line.strip().split("\t")[key_index])
                # if re.findall("(?<=_).+?(?=\])", headers[key_index]) == "ncbi":
                # if re.search("NCBI", headers[key_index]):
                #     gene.ncbi_gene_id = line.strip().split("\t")[key_index]
                #     id = True
                # if re.search("NCBI", value[0]):
                #     print(key_index, value)
                if re.search("HGNC", value[0]):
                    gene.hgnc_symbol = line.strip().split("\t")[key_index]
                    id = True
                if re.search("GenePanels_Symbol", value[0]):
                    p_symbol.symbol = line.strip().split("\t")[key_index]
                if re.search("Aliases", value[0]) and id == True:
                    aliases = []
                    # print(line.strip().split("\t")[key_index].split("|"))
                    for alias in line.strip().split("\t")[key_index].split(
                            "|"):
                        al = Alias(hgnc_symbol=alias)
                        aliases.append(al)
                    # print(line.strip().split("\t")[key_index])
                # t.append(aliases)
                value.append(gene)
            fi = FileInfo(gene=gene, alias=aliases, p_symbol=p_symbol)
            tt.append(fi)  # kan niet in 1 regel
            data_test.append(gene)

            # f = FileInfo(gene=gene)

            # t.append(FileInfo(gene=gene))

    print(len(data.get(1)))
    print(len(data_test))
    # print(data_test)
    # print(t)
    print(len(tt))
    print(tt)
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


@dataclass
class FileInfo:
    gene: Gene
    # panel: Genepanel
    # alias: List[Alias]
    p_symbol: GenepanelSymbol
    alias: list = field(default_factory=list)


def main():
    # headers = ["GeneID_NCBI", "Symbol_HGNC", "Aliases"]
    headers = ["GeneID_NCBI", "Symbol_HGNC", "Aliases", "GenePanels_Symbol"]
    # headers = ["GeneID_NCBI", "Symbol_HGNC", "Aliases", "GenePanels_Symbol"]
    file = "/Users/lean/Documenten/School/Flask/Course8_project/Parser/GenPanelOverzicht_edited.csv"
    file = "/Users/lean/Documenten/School/Flask/Course8_project/Parser/GenPanelOverzicht_DG-3.1.0_HAN_original_tsv.txt"
    reader(file, headers)
    # f = open(file)
    # print(f.read())


main()
