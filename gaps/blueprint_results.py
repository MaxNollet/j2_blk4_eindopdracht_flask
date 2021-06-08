from flask import Blueprint, render_template, url_for
from sqlalchemy import select

from gaps.models import *

blueprint_results = Blueprint("blueprint_results", __name__)

column_name_converter = {"hgnc_symbol": "HGNC symbol", "doi": "DOI", "name": "Journal"}


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
    # retrieved_values = list()
    # for rowproxy in query_results:
    #     row = dict()
    #     for field, value in zip(rowproxy.fields, rowproxy.data):
    #         row[field] = value
    #     retrieved_values.append(row)
    # for jup in retrieved_values:
    #     print(jup)
    return render_template("template_results.html")
    # search = request.args.get("search")
    # return render_template("template_query_builder.html", active="home", search=search)

# shift command r = hard refresh cache/ en zooi
