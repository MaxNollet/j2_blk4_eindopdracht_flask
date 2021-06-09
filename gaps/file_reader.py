import re


def file_reader(bestandsnaam):
    file = open(bestandsnaam, "r", encoding="utf8")
    set_genes_list = set()
    symbolen = "[, :/-]"
    for line in file:
        if not re.search(symbolen, line):
            set_genes_list.add(line.strip())
    print(set_genes_list)
    return set_genes_list
