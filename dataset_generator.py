import os


def main():
    return 0


if __name__ == "__main__":
    main()


class GenePanelParser:
    """A class which reads a gene panel file and stores
       the HGNC-symbols and aliases. These results can
       be retrieved by using corresponding getter-methods.
    """
    def __init__(self, file_name: str) -> None:
        """A method which initializes the object.

        Input = name of the file to be parsed (str).
        Output = none (None).
        """
        self.__file_name = None
        self.__symbols = set()
        self.__aliases = set()

        self.set_file_name(file_name)

    def read_file(self) -> None:
        """A method which reads a file and stores HGNC-symbols
           and aliases in separate datastructures which can be
           accessed via corresponding getter-methods.
        """
        if self.__file_name is not None:
            pass
        else:
            raise NoFileEntered
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


class IncorrectFileName(Exception):
    """Exception raised when an empty string as a filename
       was entered.
    """
    pass


class NoFileEntered(Exception):
    """Exception raised when a filename was entered as
       None-type."""
    pass
