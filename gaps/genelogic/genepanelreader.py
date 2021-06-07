from dataclasses import dataclass, field
from os import path
from re import compile, Pattern
from typing import List, Tuple, Dict


@dataclass
class GenepanelContent:
    """A class used for storing values. The stored values
       can be used to fill a database and the class can be
       used by readers/parsers to store values in a nice
       and structured way.
    """
    all_genes: List[dict] = field(default_factory=list)
    all_aliases: List[dict] = field(default_factory=list)
    all_genepanel_symbols: List[dict] = field(default_factory=list)
    all_genepanels: List[dict] = field(default_factory=list)
    all_inheritace_types: List[dict] = field(default_factory=list)
    relation_gene_alias: List[dict] = field(default_factory=list)
    relation_gene_genepanel: List[dict] = field(default_factory=list)
    relation_genepanel_inheritance: List[dict] = field(default_factory=list)


class GenepanelReader:
    """A class which has logic to read a file containing
       genepanel-data. The data read from the file will be
       made available in a GenepanelContent-object so other
       functions/methods can access values in a structured way.
    """
    _delimiter: str = "\t"
    _column_names: Tuple[str] = ("GeneID_NCBI", "Symbol_HGNC", "Aliases",
                                 "GenePanels_Symbol", "GenePanel")
    _table_names = {"alias", "gene", "genepanel", "genepanel_symbol",
                    "inheritance_type"}
    _pattern_inheritance_fix: Pattern = compile(r"(?<=[a-zA-Z0-9]);(?=[a-zA-Z0-9])")
    _pattern_find_genepanel_abbreviation: Pattern = compile(r"[a-zA-Z]+?(?=\s)")
    _pattern_find_inheritance_types: Pattern = compile(r"(?<=\()(.+?)(?=\))")
    _pattern_variable_splitter: Pattern = compile(r"\W+")

    def __init__(self, filename: str = None, auto_read: bool = True) -> None:
        """Constructor of the object, initializes a cache, saves
           a filename if given and starts automatically reading a
           file when prompted to.

        :param filename Name of the file to read (str).
        :param auto_read Start automatically reading a file or not (bool).
        """
        self._column_indexes: Dict[str, int] = dict()
        self._cached_values: Dict[str, set] = dict()
        self._data: GenepanelContent = GenepanelContent()
        self._filename: str = str()
        for table_name in self._table_names:
            self._cached_values[table_name] = set()
        if filename:
            self.set_filename(filename)
        if auto_read and self._filename is not None:
            self.read_file()

    def read_file(self, filename: str = None) -> None:
        """A method which reads a file containing data of genepanels.
           When done, all values are stored in a GenepanelContent-
           object which is available through the get_genepanel_data-
           method.

        :param filename Name of the file to read (str).
        """
        if filename:
            self.set_filename(filename)
        with open(self._filename, "r", encoding="UTF-8") as file:
            # Skip empty head of file if applicable.
            line = file.readline().strip()
            while line == "":
                line = file.readline().strip()
            self._get_column_indexes(line)
            # Extract values from the rest of the file.
            for line in file:
                stripped_line = line.strip()
                if stripped_line != "":
                    values = stripped_line.split(self._delimiter)
                    genepanel_symbol = self._extract_genepanel_symbol(values)
                    gene_symbol = self._extract_gene(values, genepanel_symbol)
                    self._extract_aliases(values, gene_symbol)
                    self._extract_genepanels(values, gene_symbol)

    def _get_column_indexes(self, line: str) -> None:
        """A method which column indexes of column names. These
           indexes are used when extracting values from desired
           columns.

        :param line First line of the file containing column names (str).
        """
        values = line.strip().split(self._delimiter)
        for column_name in self._column_names:
            try:
                self._column_indexes[column_name] = values.index(column_name)
            except ValueError:
                raise GenepanelColumnNotFound(column_name)
        return None

    def _extract_genepanel_symbol(self, values: List[str]) -> str:
        """A method which extracts the symbol which is used to display
           the gene for a certain genepanel.

        :param values One row of all columns from the file (List[str]).
        :return Symbol used for displaying the gene in a certain
                genepanel (str).
        """
        genepanel_symbol = values[self._column_indexes["GenePanels_Symbol"]].strip()
        cache = self._cached_values["genepanel_symbol"]
        if genepanel_symbol != "" and genepanel_symbol not in cache:
            cache.add(genepanel_symbol)
            self._data.all_genepanel_symbols.append({"symbol": genepanel_symbol})
        return genepanel_symbol

    def _extract_gene(self, values: List[str], genepanel_symbol: str) -> str:
        """A method which extracts all values needed for inserting
           a gene into the database.

        :param values One row of all columns from the file (List[str]).
        :param genepanel_symbol Symbol of the genepanel (str).
        :return Symbol of the gene (str).
        """
        gene_symbol = values[self._column_indexes["Symbol_HGNC"]].strip()
        cache = self._cached_values["gene"]
        if gene_symbol not in cache:
            cache.add(gene_symbol)
            gene_id = values[self._column_indexes["GeneID_NCBI"]].strip()
            self._data.all_genes.append({"ncbi_gene_id": gene_id, "hgnc_symbol": gene_symbol,
                                         "genepanel_symbol_id": genepanel_symbol,
                                         "in_genepanel": True})
        return gene_symbol

    def _extract_aliases(self, values: List[str], gene_symbol: str) -> None:
        """A method which extracts all aliases which belong to a
           specific gene.

        :param values One row of all columns from the file (List[str]).
        :param gene_symbol Symbol of the gene which these aliases
                           belong to (str).
        """
        aliases = values[self._column_indexes["Aliases"]].strip()
        cache = self._cached_values["alias"]
        if aliases != "":
            if "|" in aliases:
                for value in aliases.split("|"):
                    alias = value.strip()
                    if alias != "" and alias not in cache:
                        cache.add(alias)
                        self._data.all_aliases.append({"hgnc_symbol": alias})
                    self._data.relation_gene_alias.append({"alias_id": alias,
                                                           "gene_id": gene_symbol})
            else:
                if aliases not in cache:
                    cache.add(aliases)
                    self._data.all_aliases.append({"hgnc_symbol": aliases})
                self._data.relation_gene_alias.append({"alias_id": aliases,
                                                       "gene_id": gene_symbol})
        return None

    def _extract_genepanels(self, values: List[str], gene_symbol: str) -> None:
        """A method which extracts all genepanels with all inheritance
           types where a specific gene belongs in.

        :param values One row of all columns from the file (List[str]).
        :param gene_symbol Symbol of the gene which these genepanels
                           belong to (str).
        """
        genepanels = values[self._column_indexes["GenePanel"]].strip()
        cache_genepanel = self._cached_values["genepanel"]
        cache_inheritance = self._cached_values["inheritance_type"]
        if genepanels != "":
            line_fix = self._fix_delimeter_inconsistency(genepanels)
            for value in line_fix.split(";"):
                genepanel = self._pattern_find_genepanel_abbreviation.findall(value)[0]
                if genepanel != "" and genepanel not in cache_genepanel:
                    cache_genepanel.add(genepanel)
                    self._data.all_genepanels.append({"abbreviation": genepanel})
                inheritance_types = self._pattern_find_inheritance_types.findall(value)[0]
                if inheritance_types != "" and genepanel != "":
                    for inheritance in self._pattern_variable_splitter.split(inheritance_types):
                        if inheritance != "" and inheritance not in cache_inheritance:
                            cache_inheritance.add(inheritance)
                            self._data.all_inheritace_types.append({"type": inheritance})
                        self._data.relation_genepanel_inheritance.append({"genepanel_id": genepanel,
                                                                          "inheritance_type_id": inheritance})
                self._data.relation_gene_genepanel.append({"genepanel_id": genepanel,
                                                           "gene_id": gene_symbol})
        return None

    @classmethod
    def _fix_delimeter_inconsistency(cls, values: str):
        """A method which cleans up human-inconsistencies when specifying
           inheritance types for genepanels. All inheritance types are
           delimited by a comma-sign (,) but sometimes a semicolon (;) is
           used. This messes with the way genepanel abbreviations can be
           distinguished from one and another.

        :param values Column containing genepanel abbreviations with
                      inheritance-types to fix (str).
        :return Fixed-up values (str).
        """
        return cls._pattern_inheritance_fix.sub(",", values)

    def set_filename(self, filename: str) -> None:
        """A method which validates a filename and saves the filename
           for later uses when reading the file.

        :param filename Name of the file to read (str).
        :raises FileNotFoundError If the file can't be found.
        """
        filename = filename.strip()
        if path.exists(filename):
            self._filename = filename
        else:
            raise FileNotFoundError
        return None

    def get_genepanel_data(self) -> GenepanelContent:
        """A method used to retrieve the read values from a file with
           genepanel-information.

        :return All read values from a genepanel-file (GenepanelContent).
        """
        return self._data


class GenepanelColumnNotFound(Exception):
    """An exception when a specific column cannot be found in the
       file that had been uploaded by the user. Instead of
       providing a 'ValueError' this custom exception aims
       to help specify which column is missing from the file.
    """
    def __init__(self, column_name: str):
        super().__init__(f"The column '{column_name}' was not found in the file "
                         f"that had been uploaded! Could not update genepanels.")
