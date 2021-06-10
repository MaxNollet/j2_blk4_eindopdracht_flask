from typing import Tuple

from sqlalchemy import select

from gaps.query_builder.table_field_describer import TableFieldDescriber
from gaps.query_builder.table_join_describer import TableJoinDescriber
from gaps.models import Query


class SelectStatementBuilder(TableFieldDescriber, TableJoinDescriber):
    """A class which contains logic to build a dynamic query
       based on selected fields from the database.
    """
    def __init__(self, fields: Tuple[str, ...], query_id: str):
        """Constructor of the object, joins all tables and stores
           the generated query inside the object.

        :param fields All fields to put in the select statement (Tuple[str, ...]).
        :param query_id ID of the query which matched genes (str).
        """
        TableFieldDescriber.__init__(self)
        TableJoinDescriber.__init__(self)
        self.__statement = None
        self.__joined_tables = set()
        select_columns = self.__get_select_columns(fields)
        select_statement = self.__create_select_statement(select_columns)
        filtered_statement = self.__apply_filter(select_statement, query_id)
        self.__statement = filtered_statement

    def __get_select_columns(self, columns: Tuple[str, ...]) -> tuple:
        """A method which extracts all valid fields available to select
           and returns valid fields.

        :param columns Column/fields to select from the database (Tuple[str, ...])).
        :return All valid fields which can be selected from the database (tuple).
        """
        select_columns = list()
        for column in columns:
            column_name = self.column_name_describer.get(column)
            if column_name:
                select_columns.append(column_name)
        return tuple(select_columns)

    def __create_select_statement(self, field_names: tuple) -> select:
        """A method which processes all valid fields which have to be
           selected from the database and returns a valid query based
           on all fields.

        :param field_names Fields to put in the select statement (tuple).
        :return Select statement for all fields (Select).
        """
        fields_to_select = [element.get("column") for element in field_names]
        statement = select(fields_to_select)
        for field in field_names:
            statement = self.__join_table(field.get("table"), statement)
        return statement

    def __join_table(self, table_name: str, statement: select) -> select:
        """A method which joins required tables so all fields can
           be selected from the correct table in the database.

        :param table_name Name of the table to join (str).
        :param statement Select statement to append join-clause (Select).
        :return Select-statement with joins for required tables (Select).
        """
        if table_name not in self.__joined_tables and table_name in self.table_joins:
            self.__joined_tables.add(table_name)
            # print(table_name)  # eg gene, article_gene, disease, etc.
            join_table = self.table_joins.get(table_name)
            required_joins: list = join_table.get("required")
            if required_joins:
                for required_join in required_joins:
                    statement = self.__join_table(required_join, statement)
            table = join_table.get("table")
            id1 = join_table.get("id1")
            id2 = join_table.get("id2")
            statement = statement.join(table, onclause=id1 == id2, isouter=True)
        return statement

    def __apply_filter(self, statement: select, query_id: str):
        """A method which applies a filter in order to retrieve
           only genes which were found with a specific query.

        :param statement Statement selecting all fields (Select).
        :param query_id Query ID to filter the select statement (str).
        :return Filtered select-statement (Select).
        """
        statement = self.__join_table(Query.__tablename__, statement)
        return statement.where(self.column_name_describer.get("Query ID").get("column") == query_id)

    def __order_by(self, statement: select, order: str):
        pass

    def get_statement(self) -> select:
        """A method which returns the statement to select desired
           fields from the database.

        :return Generated select statement (Select).
        """
        return self.__statement
