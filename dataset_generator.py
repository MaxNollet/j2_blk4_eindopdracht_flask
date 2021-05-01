import os


def main():
    parser = GenePanelParser(file_name="GenPanelOverzicht_DG-3.1.0_HAN_original.tsv")
    print(len(parser.get_symbols()))
    print(len(parser.get_aliases()))
    return 0


class GenePanelParser:
    """A class which reads a gene panel file and stores
       the HGNC-symbols and aliases. These results can
       be retrieved by using corresponding getter-methods.
    """
    def __init__(self, file_name: str = None) -> None:
        """A method which initializes the object.

        Input = name of the file to be parsed (str).
        Output = none (None).
        """
        self.__file_name = None
        self.__symbols_column_name = "Symbol_HGNC"
        self.__aliases_column_name = "Aliases"
        self.__symbols = None
        self.__aliases = None
        # Set filename if given
        if file_name is not None:
            self.set_file_name(file_name)
            if self.__file_name is not None:
                self.parse_file()

    def parse_file(self) -> None:
        """A method which reads a file and stores HGNC-symbols
           and aliases in separate datastructures which can be
           accessed via corresponding getter-methods.
        """
        self.__symbols = set()
        self.__aliases = set()
        if self.__file_name is not None:
            with open(self.__file_name, "r") as file:
                line = file.readline()
                while line == "":
                    line = file.readline()
                # Empty head of file skipped, determining column indexes.
                header = line.split("\t")
                column_symbol = header.index(self.__symbols_column_name)
                column_alias = header.index(self.__aliases_column_name)
                for line in file:
                    splitted_line = line.split("\t")
                    self.__add_symbol(splitted_line[column_symbol])
                    [self.__add_alias(alias) for alias in splitted_line[column_alias].split("|")]
        else:
            raise NoFileEntered
        return None

    def __add_symbol(self, symbol: str) -> None:
        """A method which adds a symbol to the set of symbols
           and strips excessive whitespaces around the symbol.

        Input = symbol of a gene (str).
        """
        stripped_string = symbol.strip()
        if stripped_string != "":
            self.__symbols.add(stripped_string)
        return None

    def __add_alias(self, alias: str) -> None:
        """A method which adds a symbol to the set of aliases
           and strips excessive whitespaces around the alias.

        Input = alias of a gene (str).
        """
        stripped_alias = alias.strip()
        if stripped_alias != "":
            self.__aliases.add(stripped_alias)
        return None

    def set_file_name(self, file_name: str) -> None:
        """Setter for the filename to be parsed later. Performs
           some checks if the filename is not an empty string
           and if the file exists.

        Input = name of the file to be parsed (str).
        Output = none (None).
        """
        if file_name is not None:
            stripped_file_name = file_name.strip()
            if stripped_file_name != "":
                if os.path.isfile(stripped_file_name):
                    self.__file_name = stripped_file_name
                else:
                    raise FileNotFoundError
            else:
                raise IncorrectFileName
        else:
            raise NoFileEntered
        return None

    def get_symbols(self) -> tuple:
        """A method which returns a sorted tuple of all
           unique symbols parsed from the file. If no file
           has been parsed yet, an empty tuple is returned.

        Output = all unique symbols parsed form the file (tuple).
        """
        if self.__symbols is not None:
            return tuple(sorted(self.__symbols))
        else:
            return tuple()

    def get_aliases(self) -> tuple:
        """A method which returns a sorted tuple of all
           unique aliases parsed from the file. If no file
           has been parsed yet, an empty list is returned.

        Output = all unique aliases parsed from the file (tuple).
        """
        if self.__aliases is not None:
            return tuple(sorted(self.__aliases))
        else:
            return tuple()


class IncorrectFileName(Exception):
    """Exception raised when an empty string as a filename
       was entered.
    """
    pass


class NoFileEntered(Exception):
    """Exception raised when a filename was entered as
       None-type."""
    pass


if __name__ == "__main__":
    main()
