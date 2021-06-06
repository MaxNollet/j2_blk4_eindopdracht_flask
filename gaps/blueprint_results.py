from flask import Blueprint, render_template, url_for
from sqlalchemy import select

from gaps.models import *

blueprint_results = Blueprint("blueprint_results", __name__)


@blueprint_results.route("/results")
def results():
    """A function which handles requests to the '/'-
       route for the webapp.

    :return Rendered template of 'homepage.html' (str).
    """
    joins = select([Gene.hgnc_symbol, Article.doi, Journal.name])\
        .join(t_article_gene, Gene.id == t_article_gene.c.gene_id)\
        .join(Article, Article.id == t_article_gene.c.article_id)\
        .join(Journal)\
        .join(t_query_gene, Gene.id == t_query_gene.c.gene_id)\
        .join(Query, Query.id == t_query_gene.c.query_id)\
        .where(Query.id == "df56db6d-ce3a-4aa3-b91e-48c455355f95")\
        .order_by(Gene.hgnc_symbol)
    print(joins)
    query_results = db.session.execute(joins)
    retrieved_values = list()
    for rowproxy in query_results:
        row = dict()
        for field, value in zip(rowproxy._fields, rowproxy._data):
            row[field] = value
        retrieved_values.append(row)
    for jup in retrieved_values:
        print(jup)
    # for row in results:
    #     print(row.__dict__)
    return render_template("template_results.html")
    # search = request.args.get("search")
    # return render_template("template_query_builder.html", active="home", search=search)

# shift command r = hard refresh cache/ en zooi
