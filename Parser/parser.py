from dataclasses import dataclass, field
from gaps.models import Gene, Alias, Genepanel, GenepanelSymbol, \
    InheritanceType
import re


def reader(file, headers):
    """
    Reads the GenePanel file, gets the headers from the first row
    :param file: GenPanelOverzicht | csv file | first row needs to contain
    the names of the columns
    :param headers: names of the columns to get the data from
    :return:
    """
    print("reading file")
    data = {}  # index and name column
    file_list = []  # Contains objects from FileInfo
    with open(file, mode="r", encoding="utf-8-sig") as f:
        # https://stackoverflow.com/questions/17912307/
        # u-ufeff-in-python-string
        # data = [i for i, j in enumerate(f.readline().strip()
        # .split(";")) if "GenePanels" in j]
        for i, j in enumerate(f.readline().strip().split("\t")):
            # i is the index and j is the value
            for h in headers:
                if h == j:  # headers needs to match with column
                    # data.update({i: j})  # hier staat er \ufeff voor
                    data[i] = [j]
        for line in f:
            gene = Gene(in_genepanel=True)
            p_symbol = GenepanelSymbol()
            for key_index, value in data.items():  # beter naam voor
                if re.search("HGNC", value[0]):
                    gene.hgnc_symbol = line.strip().split("\t")[key_index]
                if re.search("GenePanels_Symbol", value[0]):
                    p_symbol.symbol = line.strip().split("\t")[key_index]
                if re.search("Aliases", value[0]):
                    aliases = []
                    # print(line.strip().split("\t")[key_index].split("|"))
                    for alias in line.strip().split("\t")[key_index].split(
                            "|"):
                        al = Alias(hgnc_symbol=alias)
                        aliases.append(al)
                if "GenePanel" == value[0]:
                    panels = []  # list with 'LENGTE','MR', 'SCHISIS' eg
                    # currently not in use
                    ihh_list = []  # list with inheritance
                    combi_panel = []  # combi [OMIM, AR, AD] example
                    # print(line.strip().split("\t")[key_index]) # eg
                    haken = re.findall(f"(?<=\().+?(?=\))",
                                       line.strip().split("\t")[key_index])
                    print(haken, "HAKEN")
                    if re.search(";", str(haken)):
                        for h in haken:
                            print(h, " h")
                            testje = h.replace(";", ",")
                            # print(line.strip().split("\t")[key_index])  # eg
                            test = re.sub(h, testje,
                                          line.strip().split("\t")[key_index])
                            test = test.split(";")
                            print(test, "testen")
                            # print(test)
                    else:
                        # print(line.strip().split("\t")[key_index])  # eg
                        test = line.strip().split("\t")[key_index].split(";")
                    for te in test:
                        een_genpanel = []
                        p = re.findall(f"(?<=\().+?(?=\))",
                                       te)
                        print(test)
                        print(p, "pp")
                        # print(p)
                        # the inheritance between ()
                        # print(p)
                        # test = line.strip().split("\t")[key_index].split(";")
                        # print(test)
                        print(te, "te")
                        if len(p) >= 1:
                            # print(p, "P")
                            for ih in p:  # p = ['AD', 'AD',
                                # 'AD;UK,AR,AD,XL', 'AD']
                                # print(ih, " IH")
                                if ";" or "," in str(ih):  # ih 'AD;UK,XL'
                                    ihh = str(ih).replace(";", ",").split(",")
                                    for j in ihh:
                                        print(j, "jjjjj")
                                        een_genpanel.append(
                                            InheritanceType(type=j))
                                        # ihh_list.append(InheritanceType(type=j))
                                else:
                                    for j in ih:  # example 'AD'
                                        een_genpanel.append(
                                            InheritanceType(type=j))
                                        # ihh_list.append(InheritanceType(type=j))

                        t = re.sub(f"(?<=\().+?(?=\))", "", te)
                        # t = re.sub(f"(?<=\().+?(?=\))", "", str(line.strip().split("\t")[key_index]))
                        k = t.replace(" ()", "").replace("\"", "").split(";")
                        # k = t.replace(" ()", "").replace("\"", "").split(";")
                        for a in k:  # k contains ['HEMOS', 'OMIM'] example
                            gp = Genepanel(abbreviation=a)
                            panels.append(gp)
                            een_genpanel.append(gp)
                        combi_panel.append(een_genpanel)
                        # print(re.findall(f".+?(?=\s\()", t))
            # fi = FileInfo(gene=gene, alias=aliases, p_symbol=p_symbol,
            #               panel=panels, p_inheritance=ihh_list)
            fi = FileInfo(gene=gene, alias=aliases, p_symbol=p_symbol,
                          panel=combi_panel)
            file_list.append(fi)  # kan niet in 1 regel
    print("complete")
    # print(file_list)
    print(file_list[len(file_list) - 6])
    print(file_list[len(file_list) - 1])
    return tuple(file_list)


@dataclass
class FileInfo:
    gene: Gene
    p_symbol: GenepanelSymbol
    panel: list = field(default_factory=list)  # includes AR, AD en afkorting
    # p_inheritance: list = field(default_factory=list)
    alias: list = field(default_factory=list)


def main():
    # headers = ["GeneID_NCBI", "Symbol_HGNC", "Aliases"]
    headers = ["GeneID_NCBI", "Symbol_HGNC", "Aliases", "GenePanels_Symbol",
               "GenePanel"]
    file = "/Users/lean/Documenten/School/Flask/Course8_project/" \
           "Parser/GenPanelOverzicht_DG-3.1.0_HAN_original_tsv.txt"
    try:
        reader(file, headers)
    except FileNotFoundError:
        print("File not found! ->", file)


main()
