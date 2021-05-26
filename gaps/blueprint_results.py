from flask import Blueprint, render_template, url_for

blueprint_results = Blueprint("blueprint_results", __name__)


@blueprint_results.route("/results")
def results():
    """A function which handles requests to the '/'-
       route for the webapp.

    :return Rendered template of 'homepage.html' (str).
    """
    return render_template("template_results.html")
    # search = request.args.get("search")
    # return render_template("template_query_builder.html", active="home", search=search)

# shift command r = hard refresh cache/ en zooi
