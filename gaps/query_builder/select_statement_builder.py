from typing import Tuple

from sqlalchemy import select

from gaps.query_builder.table_field_describer import TableFieldDescriber
from gaps.query_builder.table_join_describer import TableJoinDescriber
from gaps.models import Query


class SelectStatementBuilder(TableFieldDescriber, TableJoinDescriber):
    def __init__(self, fields: Tuple[str, ...], query_id: str):
        TableFieldDescriber.__init__(self)
        TableJoinDescriber.__init__(self)
        self.__statement = None
        self.__joined_tables = set()
        select_columns = self.__get_select_columns(fields)
        select_statement = self.__create_select_statement(select_columns)
        filtered_statement = self.__apply_filter(select_statement, query_id)
        self.__statement = filtered_statement

    def __get_select_columns(self, columns: Tuple[str, ...]) -> tuple:
        select_columns = list()
        for column in columns:
            column_name = self.column_name_describer.get(column)
            if column_name:
                select_columns.append(column_name)
        return tuple(select_columns)

    def __create_select_statement(self, field_names: tuple) -> select:
        fields_to_select = [element.get("column") for element in field_names]
        statement = select(fields_to_select)
        for field in field_names:
            statement = self.__join_table(field.get("table"), statement)
        return statement

    def __join_table(self, table_name: str, statement: select) -> select:
        if table_name not in self.__joined_tables and table_name in self.table_joins:
            self.__joined_tables.add(table_name)
            print(table_name)
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
        statement = self.__join_table(Query.__tablename__, statement)
        # return statement.where(Query.id == query_id)
        return statement.where(self.column_name_describer.get("Query ID").get("column") == query_id)

    def order_by(self, statement: select, order: str):
        pass

    def get_statement(self) -> select:
        return self.__statement
