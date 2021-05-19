from flask import Blueprint, render_template, request

blueprint_home = Blueprint("blueprint_home", __name__)


@blueprint_home.route("/")
def homepage():
    """A function which handles requests to the '/'-
       route for the webapp.

    :return Rendered template of 'homepage.html' (str).
    """
    search = request.args.get("search")
    return render_template("template_query.html", active="home", search=search)
# shift command r = hard refresh cache/ en zooi
