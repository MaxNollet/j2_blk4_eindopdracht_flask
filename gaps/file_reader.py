import re

def main():
    bestandsnaam = "file_test.txt"
    file_reader(bestandsnaam)

def file_reader(bestandsnaam):
    file = open(bestandsnaam, "r", encoding="utf8")
    genesymbol_list = []
    symbolen = "[, :/-]"
    for line in file:
        if not re.search(symbolen, line):
            genesymbol_list.append(line.strip())
    print(genesymbol_list)

main()
