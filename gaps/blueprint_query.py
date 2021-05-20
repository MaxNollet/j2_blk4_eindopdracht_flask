from flask import Blueprint, render_template

from gaps.genelogic import DatabaseInserter

blueprint_query = Blueprint("blueprint_query", __name__)


@blueprint_query.route("/query")
def query():
    """A function which handles requests of the '.query'-
       route for the webapp.

    :return Rendered template of 'query.html' (str).
    """
    DatabaseInserter.updateGenpanel()
    return render_template("query.html", active="query_input")

# shift command r = hard refresh cache/ en zooi
