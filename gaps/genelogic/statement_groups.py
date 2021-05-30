from typing import Dict

from sqlalchemy.dialects.postgresql import Insert
from sqlalchemy.sql import Select


class StatementGroups:
    """A class which contains logic for adding statements
       to a central mapper. Statements are categorized per
       table name and are classified by type (insert, select,
       delete, etc.). This mapper can be extended by other
       classes so that other classes can retrieve statements
       for a desired table.
    """
    _statement_groups: Dict[str, Dict[str, object]] = {}

    @staticmethod
    def add_to_group(statement, table: str, column_as_key: str = None) -> None:
        """A method which stores statements inside the central
           mapper and categorizes them based on the type of the
           statement.

        :param statement Statement to store in the central mapper.
        :param table Name of the table the statement belongs to (str).
        :param column_as_key Name of the column to select values from
                             to be used as keys for dictionaries (str).
        """
        if StatementGroups._statement_groups.get(table) is None:
            StatementGroups._statement_groups[table] = dict()
        if column_as_key is not None:
            StatementGroups._statement_groups[table]["column_as_key"] = column_as_key
        if isinstance(statement, Insert):
            StatementGroups._statement_groups[table]["insert"] = statement
        elif isinstance(statement, Select):
            StatementGroups._statement_groups[table]["select"] = statement


def statement_group(table: str, column_as_key: str = None):
    """A decorator which adds statements to a central mapper
       where statements can be retrieved from in an orderly
       manner.

    :param table Name of the table to which the statement belongs to (str).
    :param column_as_key Name of the column to select values from to be
                         used as a key for dictionaries (str).
    :return Decorated function of the statement (Function).
    """

    def decorator_add(statement_function):
        """A function which extracts the statement from a function
           and stores the statement in a central mapper.

        :param statement_function Function containing the statement
                                  to be stored in the mapper (Function).
        :return Original function containing the statement (Function).
        """
        StatementGroups.add_to_group(statement_function(), table, column_as_key)
        return statement_function

    return decorator_add
