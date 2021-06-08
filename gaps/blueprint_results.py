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
    joins = select([Gene.hgnc_symbol, Article.doi, Article.journal_id])\
        .join(t_article_gene, Gene.id == t_article_gene.c.gene_id)\
        .join(Article, Article.id == t_article_gene.c.article_id)\
        .join(Journal)\
        .join(t_query_gene, Gene.id == t_query_gene.c.gene_id)\
        .join(Query, Query.id == t_query_gene.c.query_id)\
        .where(Query.id == query_id)\
        .order_by(Gene.hgnc_symbol)
    print(joins)

    query_results = db.session.execute(joins)
    for row in query_results:
        print(row)
    return render_template("template_results.html")


class SelectStatementBuilder:
    column_name_getter = {"HGNC symbol": Gene.hgnc_symbol, "DOI": Article.doi, "Journal": Article.journal_id}
    column_name_converter = {"hgnc_symbol": "HGNC symbol", "doi": "DOI", "name": "Journal"}

    def select_statement_builder(self, columns: tuple):
        for column in columns:
            pass
        return None
