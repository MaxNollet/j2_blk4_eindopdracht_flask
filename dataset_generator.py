import os
import random
import matplotlib.pyplot as plot


def main():
    symbols = SymbolParser(file_name="hgnc_complete_set.txt",
                           symbols_column_name="symbol",
                           aliases_column_name="alias_symbol")
    # GraphGenerator((symbols.get_aliases(), symbols.get_symbols()),
    #                ("Aliases", "Symbols"),
    #                "Lengths of aliases and symbols from HGNC")
    words = WordListReader(file_name="words.txt")
    print(f"Symbols: {len(symbols.get_symbols())}")
    print(f"Aliases: {len(symbols.get_aliases())}")
    print(f"Words: {len(words.get_words())}")
    equalizer = DatasetEqualizer(symbols=symbols, words=words)
    return 0


def pretty_print(categorized: dict):
    # print(sorted(categorized.keys()))
    for key in sorted(categorized.keys()):
        print(f"{key}: {len(categorized.get(key))}")
    return None


class SymbolParser:
    """A class which reads a gene panel file and stores
       the HGNC-symbols and aliases. These results can
       be retrieved by using corresponding getter-methods.
    """

    def __init__(self, file_name: str, auto_parse: bool = True,
                 symbols_column_name: str = "Symbol_HGNC",
                 aliases_column_name: str = "Aliases") -> None:
        """A method which initializes the object.

        Input = -name of the file to be parsed (str).
                -indication to start parsing automatically (bool).
                -column name where symbols are stored (sr).
                -column name where aliases are stored (str).
        Output = none (None).
        """
        self.__file_name = None
        self.__symbols_column_name = symbols_column_name
        self.__aliases_column_name = aliases_column_name
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
            with open(self.__file_name, "r", encoding="UTF-8") as file:
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


class WordListReader:
    """A class which generates a list of words and filters
       this list according to requirements.
    """

    def __init__(self, file_name: str, auto_parse: bool = True) -> None:
        """A method which constructs the object and starts
           to generate a list of words automatically if the
           variable 'auto_generate' is set to True.

        Input = indication to start generating a list of words automatically (bool).
        """
        self.__filename = file_name
        self.__words = None
        # start generating a list of words automatically.
        if auto_parse:
            self.__words = self.read_words()

    def read_words(self) -> set:
        """A method which reads all words from a file and stores
           all unique words in a set.

        Output = all unique words from the specified file (set).
        """
        if self.__filename is not None and self.__filename != "":
            words = set()
            with open(self.__filename, "r", encoding="UTF-8") as file:
                for line in file:
                    word = line.strip()
                    if word != "":
                        words.add(word)
            return words
        else:
            raise NoFileEntered

    def get_words(self) -> set:
        """a method which returns a list of all unique
           words from the specified file.

        Output = list of all unique retrieved words (set).
        """
        if self.__words is not None:
            return self.__words
        else:
            return set()


class DatasetEqualizer:
    def __init__(self, symbols: SymbolParser,
                 words: WordListReader,
                 auto_process: bool = True,
                 equalize_dataset: bool = True,
                 output_file: str = "dataset") -> None:
        self.__symbol_parser = symbols
        self.__word_reader = words
        self.__output_file = output_file
        if auto_process:
            self.__extract_unique_variables()
            if equalize_dataset:
                self.equalize_dataset()
            # self.write_dataset()

    def __extract_unique_variables(self) -> None:
        """A method which extracts all unique symbols, aliases
           and words from the parsers and stores them inside
           the object.
        """
        symbols = self.__symbol_parser.get_symbols()
        aliases = self.__symbol_parser.get_aliases()
        words = self.__word_reader.get_words()
        self.__symbols = symbols
        self.__aliases = aliases.difference(symbols)
        self.__words = words.difference(aliases).difference(symbols)
        return None

    def equalize_dataset(self):
        """A method which uses other methods to equalize the
           provided symbols and words.
        """
        self.__determine_dataset_size()
        lengths_smaller_dataset = self.__categorize_lengths(self.__smaller_dataset)
        lengths_bigger_dataset = self.__categorize_lengths(self.__bigger_dataset)
        equalized_bigger_dataset = self.__equalize_lengths(lengths_smaller_dataset, lengths_bigger_dataset)
        return 0

    def __determine_dataset_size(self) -> None:
        """A method which compares the two datasets and categorizes
           the bigger and the smaller datasets. Based on this, the
           bigger dataset will be sampled to match the size of the
           smaller dataset.
        """
        self.__smaller_dataset = set()
        self.__bigger_dataset = set()
        size_symbols = len(self.__symbols) + len(self.__aliases)
        size_words = len(self.__words)
        if size_symbols >= size_words:
            self.__smaller_dataset.update(self.__words)
            self.__bigger_dataset.update(self.__symbols, self.__aliases)
        else:
            self.__smaller_dataset.update(self.__symbols, self.__aliases)
            self.__bigger_dataset.update(self.__words)
        return None

    def __categorize_lengths(self, dataset: set) -> dict:
        """A method which sorts all elements from a given
           dataset on the length of an element.

        Input = dataset to be sorted (set).
        Output = sorted dataset in length (dict).
        """
        lengths = dict()
        for element in dataset:
            length = len(element)
            if length not in lengths:
                lengths[length] = 1
            else:
                lengths[length] += 1
        return lengths

    def __equalize_lengths(self, smaller_dataset: dict, bigger_dataset: dict) -> dict:
        """A method which equalizes the lengths categories from
           the smaller and the bigger dataset.

        Input = -lengths of the smaller dataset (dict).
        Input = -lengths of the bigger dataset (dict).
        Output = equalized bigger dataset (dict).
        """
        to_be_remove = set()
        for element in bigger_dataset.keys():
            if element not in smaller_dataset:
                to_be_remove.add(element)
        if len(to_be_remove) > 0:
            for remove in to_be_remove:
                bigger_dataset.pop(remove)
        return bigger_dataset

    def write_dataset(self) -> None:
        """A method which writes the dataset to a file. Each line
           in the file contains one symbol or word with a
           corresponding label.
        """
        if self.__symbols is not None and self.__symbols is not None and self.__words is not None:
            with open(self.__output_file, "w", encoding="UTF-8") as file:
                for symbol in sorted(self.__symbols):
                    file.write(f"{symbol}\tsymbol\n")
                for alias in sorted(self.__aliases):
                    file.write(f"{alias}\tsymbol\n")
                for word in sorted(self.__words):
                    file.write(f"{word}\tword\n")
        return None

    @staticmethod
    def __sort_dictionary(dictionary: dict) -> list:
        return sorted(dictionary.items(), key=lambda element: len(element[1]))

    @staticmethod
    def __convert_to_lengths(values: set) -> list:
        frequencies = list()
        [frequencies.append(len(element)) for element in values]
        return sorted(frequencies)


class GraphGenerator:
    def __init__(self, dataset: tuple, labels: tuple, title: str) -> None:
        self.__dataset = dataset
        self.__types = labels
        self.__title = title
        # Variables from analysis.
        self.__lengths = None
        self.__bins = None
        self.__labels = None
        self.__numbers = None
        # Calling methods to prep the data for the graph.
        self.determine_lengths()
        self.determine_bins()
        self.numbers()
        self.histogram()

    def determine_lengths(self):
        self.__lengths = list()
        for values in self.__dataset:
            lengths = list()
            [lengths.append(len(length)) for length in values]
            self.__lengths.append(lengths)
        return None

    def determine_bins(self) -> None:
        bins = set()
        for lengths in self.__lengths:
            bins.update(lengths)
        self.__labels = sorted(bins)
        bins.add(max(bins) + 1)
        self.__bins = sorted(bins)
        return None

    def numbers(self):
        numbers = 0
        no_numbers = 0
        for dataset in self.__dataset:
            for element in dataset:
                if any(character.isdigit() for character in element):
                    numbers += 1
                else:
                    no_numbers += 1
        if numbers != 0 and no_numbers != 0:
            self.__numbers = (numbers, no_numbers)
        return None

    def histogram(self) -> None:
        plot.hist(x=self.__lengths, bins=self.__bins, stacked=True,
                  align="left", linewidth=0.5, edgecolor="black")
        plot.xticks(self.__labels)
        plot.title(self.__title)
        plot.xlabel("Length")
        plot.ylabel("Frequency")
        plot.legend(self.__types)
        # plot.tight_layout()
        # plot.savefig(f"{self.__title}.svg")
        plot.show()
        plot.close()
        return None

    def donut(self) -> None:
        plot.pie(self.__numbers, autopct="%.0f%%")
        plot.title("test")
        plot.legend(("test1", "test2"))
        return None


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


if __name__ == "__main__":
    main()
