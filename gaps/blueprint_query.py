from flask import Blueprint, render_template

blueprint_query = Blueprint("blueprint_query", __name__)


@blueprint_query.route("/query")
def homepage():
    """A function which handles requests of the '.query'-
       route for the webapp.

    :return Rendered template of 'query.html' (str).
    """
    return render_template("query.html", active="query_input")

# shift command r = hard refresh cache/ en zooi
