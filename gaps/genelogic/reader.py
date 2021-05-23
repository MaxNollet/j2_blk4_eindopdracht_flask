from dataclasses import dataclass, field
from gaps.models import Gene, Alias, Genepanel, GenepanelSymbol, \
    InheritanceType
import re
import os


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
                if re.search("GeneID_NCBI", value[0]):
                    gene.ncbi_gene_id = line.strip().split("\t")[key_index]
                if re.search("HGNC", value[0]):
                    gene.hgnc_symbol = line.strip().split("\t")[key_index]
                if re.search("GenePanels_Symbol", value[0]):
                    p_symbol.symbol = line.strip().split("\t")[key_index]
                if re.search("Aliases", value[0]):
                    aliases = []
                    for alias in line.strip().split("\t")[key_index].split(
                            "|"):
                        al = Alias(hgnc_symbol=alias, genes=[gene])
                        # al = Alias(hgnc_symbol=alias) oud
                        aliases.append(al)
                if "GenePanel" == value[0]:
                    combi_panel = []  # combi [OMIM],[ AR, AD]] example
                    haken = re.findall(f"(?<=\().+?(?=\))",
                                       line.strip().split("\t")[key_index])
                    if re.search(";", str(haken)):
                        # checks if ; in (ab; ar, xl) etc
                        for h in haken:  # loop over alle gevonden ()
                            if ";" in h:
                                repl = h.replace(";", ",")
                                # replace de ; binnen de ()
                                line_fix = re.sub(h, repl, line.strip().
                                                  split("\t")[key_index])
                                # vervang de () met ; door een ,
                                line_fix = line_fix.split(";")
                    else:
                        # print(line.strip().split("\t")x[key_index])  # eg
                        line_fix = line.strip().split("\t")[key_index].split(
                            ";")
                    for li in line_fix:  # line_fix is list from split
                        een_genpanel = []  # new for each new combi
                        all_ih = re.findall(f"(?<=\().+?(?=\))", li)
                        # the inheritance between ()
                        if len(all_ih) >= 0:  # vgm niet nodig
                            for ih in all_ih:  # p = ['AD', 'AD'] ih= UK, AR
                                # 'AD;UK,AR,AD,XL', 'AD']
                                if ";" or "," in str(ih):  # ih 'AD;UK,XL'
                                    ihh = str(ih).replace(";", ",").split(",")
                                    for j in ihh:
                                        een_genpanel.append(
                                            InheritanceType(type=j))
                                else:
                                    for j in ih:  # example 'AD'
                                        een_genpanel.append(
                                            InheritanceType(type=j))
                        t = re.sub(f"(?<=\().+?(?=\))", "", li)
                        k = t.replace(" ()", "").replace("\"", "").split(";")
                        # k = ['OMIM'] for example
                        # een_genpanel.append(Genepanel(abbreviation=k[0]))
                        test = Genepanel(abbreviation=k[0],
                                         inheritance_types=een_genpanel)
                        # combi_panel.append(een_genpanel)
                        combi_panel.append(test)
                        gene.genepanels = combi_panel
            fi = FileInfo(gene=gene, alias=aliases, p_symbol=p_symbol,
                          panel=combi_panel)
            file_list.append(fi)  # kan niet in 1 regel
    print("complete")
    print(len(file_list))
    # print(file_list[len(file_list) - 6])  # 79755 GeneID_NCBI
    # print(file_list[len(file_list) - 1])  # examples
    print(file_list[1569])  # 1571 ncbi_geneID = 8139
    return list(file_list)


@dataclass
class FileInfo:
    gene: Gene
    p_symbol: GenepanelSymbol
    panel: list = field(default_factory=list)  # includes AR, AD en afkorting
    alias: list = field(default_factory=list)


def get_reader(file):
    try:
        if os.path.exists(file):

            headers = ["GeneID_NCBI", "Symbol_HGNC", "Aliases",
                       "GenePanels_Symbol",
                       "GenePanel"]
            # file = "/Users/lean/Documenten/School/Flask/Course8_project/gaps/genelogic/GenPanelOverzicht_DG-3.1.0_HAN_original_tsv.txt"
            return reader(file, headers)
        else:
            print("Kan bestand niet vinden")
    except FileNotFoundError:
        print("File not found! ->", file)

# def main():
#     # headers = ["GeneID_NCBI", "Symbol_HGNC", "Aliases"]
#     headers = ["GeneID_NCBI", "Symbol_HGNC", "Aliases", "GenePanels_Symbol",
#                "GenePanel"]
#     file = "/gaps/genelogic/GenPanelOverzicht_DG-3.1.0_HAN_original_tsv.txt"
#     try:
#         reader(file, headers)
#     except FileNotFoundError:
#         print("File not found! ->", file)
#
#
# main()
