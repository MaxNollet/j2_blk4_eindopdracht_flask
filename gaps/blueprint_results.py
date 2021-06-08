from typing import Tuple

from flask import Blueprint, render_template, url_for
from sqlalchemy import select

from gaps.models import *

blueprint_results = Blueprint("blueprint_results", __name__)


@blueprint_results.route("/results/<query_id>")
def results(query_id):
    """A function which handles requests to the '/'-
       route for the webapp.

    :return Rendered template of 'homepage.html' (str).
    """
    joins = select([Gene.hgnc_symbol, Disease.disease])\
        .join()\
        .order_by(Gene.hgnc_symbol)
    query_results = db.session.execute(joins)
    for row in query_results:
        print(row)
    print(joins)
    # SelectStatementBuilder()
    return render_template("template_results.html")


class SelectStatementBuilder:
    column_name_getter = {"HGNC symbol": Gene.hgnc_symbol, "DOI": Article.doi, "Journal": Article.journal_id}
    column_name_converter = {"hgnc_symbol": "HGNC symbol", "doi": "DOI", "name": "Journal"}
    table_joins = {""}

    def __init__(self):
        self.statement = None
        self.select_statement_builder(("HGNC symbol", "DOI", "Journal"))
        print(self.statement)
        db.session.execute(self.statement)

    def select_statement_builder(self, columns: Tuple[str, ...]):
        select_columns = self.get_select_columns(columns)

    def get_select_columns(self, columns: Tuple[str, ...]) -> tuple:
        select_columns = list()
        for column in columns:
            column_name = self.column_name_getter.get(column)
            if column_name:
                select_columns.append(column_name)
        self.statement = select(select_columns)
        return tuple(select_columns)

    def join_tables(self):
        pass
