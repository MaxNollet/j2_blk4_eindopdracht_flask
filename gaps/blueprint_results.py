from flask import Blueprint, render_template, jsonify

from gaps.models import *
from query_builder import SelectStatementBuilder

blueprint_results = Blueprint("blueprint_results", __name__)


@blueprint_results.route("/results/<query_id>")
def results(query_id: str):
    """A function which handles requests to the '/'-
       route for the webapp.

    :return Rendered template of 'homepage.html' (str).
    """
    # Possible fields to select:
    # Alias ID, Alias Symbol, Article DOI, Article ID, Article PubMed ID,
    # Article abstract, Article publication date, Article title, Disease,
    # Disease ID, Disease MESH ID, Gene ID, Gene symbol, Genepanel ID,
    # Genepanel abbreviation, Genepanel symbol, Genepanel symbol ID,
    # In genepanel, Inheritance type, Inheritance type ID, Journal ID,
    # Journal name, NCBI gene ID, Query, Query ID.
    # See builder.get_field_names().
    fields = ("Gene symbol", "NCBI gene ID", "In genepanel", "Article title", "Article PubMed ID",
              "Article publication date", "Disease", "Disease MESH ID")
    builder = SelectStatementBuilder(fields, query_id)
    query_results = db.session.execute(builder.get_statement())

    test_results = list()
    group_by = "Gene symbol"
    group_index = fields.index(group_by)
    cache_element = None
    cache_group = None
    # for result in query_results:
    #     element = result[group_index]
    #     if element != cache_element:
    #         if cache_group is not None:
    #             print(cache_group)
    #             test_results.append(cache_group)
    #         cache_group = list()
    #     cache_element = result[group_index]
    #     for index, column in enumerate(result):
    #         if len(cache_group) <= index:
    #             cache_group.append(list())
    #         group = cache_group[index]
    #         if column not in group:
    #             group.append(column)
        # print(result)

    return render_template("template_results.html", fields=fields, results=query_results)
