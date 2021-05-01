import os
import random
import requests


def main():
    parser = GenePanelParser(file_name="GenPanelOverzicht_DG-3.1.0_HAN_original.tsv")
    print(f"Count of all symbols: {len(parser.get_symbols())}")
    print(f"Count of all aliases: {len(parser.get_aliases())}")
    print(f"Sum of symbols and aliases: {len(parser.get_symbols()) + len(parser.get_aliases())}")
    generator = WordListGenerator()
    print(f"Count of all retrieved words: {len(generator.get_word_list())}")
    DatasetGenerator(symbols=parser.get_symbols(), aliases=parser.get_aliases(), words=generator.get_word_list())
    return 0


class GenePanelParser:
    """A class which reads a gene panel file and stores
       the HGNC-symbols and aliases. These results can
       be retrieved by using corresponding getter-methods.
    """

    def __init__(self, file_name: str = None, auto_parse: bool = True) -> None:
        """A method which initializes the object.

        Input = -name of the file to be parsed (str).
                -indication to start parsing automatically (bool).
        Output = none (None).
        """
        self.__file_name = None
        self.__symbols_column_name = "Symbol_HGNC"
        self.__aliases_column_name = "Aliases"
        self.__symbols = None
        self.__aliases = None
        # Set filename if given and start parsing automatically.
        if file_name is not None:
            self.set_file_name(file_name)
            if auto_parse and self.__file_name is not None:
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

    def get_symbols(self) -> set:
        """A method which returns a set of all unique
           symbols parsed from the file. If no file has
           been parsed yet, an empty set is returned.

        Output = all unique symbols parsed form the file (set).
        """
        if self.__symbols is not None:
            return self.__symbols
        else:
            return set()

    def get_aliases(self) -> set:
        """A method which returns a set of all unique
           aliases parsed from the file. If no file has
           been parsed yet, an empty set is returned.

        Output = all unique aliases parsed from the file (set).
        """
        if self.__aliases is not None:
            return self.__aliases
        else:
            return set()


class WordListGenerator:
    """A class which generates a list of words and filters
       this list according to requirements.
    """
    def __init__(self, auto_generate: bool = True) -> None:
        """A method which constructs the object and starts
           to generate a list of words automatically if the
           variable 'auto_generate' is set to True.

        Input = indication to start generating a list of words automatically (bool).
        """
        self.__url_word_list = "http://www.mit.edu/~ecprice/wordlist.10000"
        self.__retrieved_words = None
        # start generating a list of words automatically.
        if auto_generate:
            self.generate_list()

    def generate_list(self) -> None:
        """A method which calls other methods to retrieve, filter
           and save a list of generated words. The list will be
           available through a corresponding getter.
        """
        self.__retrieved_words = self.__retrieve_words(self.__url_word_list)
        return None

    @staticmethod
    def __retrieve_words(url_word_list: str) -> set:
        """A method which retrieved a list of English words
           and saves all unique words inside the object.

        Input = url where a list of words can be downloaded (str).
        Output = unique words downloaded form the url (set).
        """
        r = requests.get(url_word_list)
        if r.status_code == 200:
            word_list = set()
            for word in str(r.text).split("\n"):
                stripped_word = word.strip()
                if stripped_word != "":
                    word_list.add(stripped_word)
            return word_list
        else:
            raise ErrorWhileDownloading(r.status_code)

    def get_word_list(self) -> set:
        """a method which returns a list of all unique
           retrieved words.

        Output = list of all unique retrieved words (set).
        """
        if self.__retrieved_words is not None:
            return self.__retrieved_words
        else:
            return set()


class DatasetGenerator:
    def __init__(self, symbols: set, aliases: set, words: set,
                 auto_process: bool = True,
                 equal_dataset: bool = True,
                 capitalize_words_half: bool = True) -> None:
        self.__symbols = symbols
        self.__aliases = aliases.difference(symbols)
        self.__words = words.difference(symbols).difference(aliases)
        if auto_process:
            self.equalize_dataset()
            pass

    def equalize_dataset(self) -> None:
        count_symbols = len(self.__symbols)
        count_aliases = len(self.__aliases)
        count_words = len(self.__words)
        if count_symbols + count_aliases > count_words:
            test = self.__equalize_lengths(self.__aliases, 5000)
            print(len(test))
        return None

    @staticmethod
    def __equalize_lengths(unequal: set, size: int) -> set:
        categorized = DatasetGenerator.__categorize_lengths(unequal)
        groups_count = len(categorized)
        group_size = round(size/groups_count)
        equalized = set()
        for _, values in DatasetGenerator.__sort_dictionary(categorized):
            count = len(values)
            if count == group_size:
                equalized.update(values)
            elif count > group_size:
                equalized.update((random.sample(values, group_size)))
            else:
                equalized.update(values)
            groups_count -= 1
            if groups_count > 0:
                group_size = round((size-len(equalized))/groups_count)
        return equalized

    @staticmethod
    def __categorize_lengths(uncategorized: set) -> dict:
        categorized = dict()
        for item in uncategorized:
            length = len(item)
            if length in categorized:
                categorized[length].add(item)
            else:
                categorized[length] = {item}
        return categorized

    @staticmethod
    def __sort_dictionary(dictionary: dict) -> list:
        return sorted(dictionary.items(), key=lambda element: len(element[1]))


class IncorrectFileName(Exception):
    """Exception raised when an empty string as a filename
       was entered.
    """
    pass


class NoFileEntered(Exception):
    """Exception raised when a filename was entered as
       None-type.
    """
    pass


class ErrorWhileDownloading(Exception):
    """Exception raised when the download of the word list
       did not return a 200 (=OK)-statuscode.
    """
    def __init__(self, status_code: int) -> None:
        """A method which constructs the object and send a
           message to the superclass to be displayed.

        Input = status code of the download (int).
        """
        self.message = f"Download exited with status code {status_code}"
        super().__init__(self.message)


if __name__ == "__main__":
    main()
