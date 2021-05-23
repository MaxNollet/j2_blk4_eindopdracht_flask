from flask import Blueprint, render_template

blueprint_query_builder = Blueprint("blueprint_query_builder", __name__)


@blueprint_query_builder.route("/query_builder")
def query_builder():
    """A function which handles requests of the '.query_builder'-
       route for the webapp.

    :return Rendered template of 'query.html' (str).
    """
    return render_template("template_query_builder.html", active_nav="query")

# shift command r = hard refresh cache/ en zooi
