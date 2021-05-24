from flask import Blueprint, redirect, url_for

blueprint_homepage = Blueprint("blueprint_homepage", __name__)


@blueprint_homepage.route("/")
def homepage():
    """A function which handles requests to the '/'-
       route for the webapp.

    :return Rendered template of 'homepage.html' (str).
    """
    return redirect(url_for("blueprint_query_builder.query_builder"))
    # search = request.args.get("search")
    # return render_template("template_query_builder.html", active="home", search=search)

# shift command r = hard refresh cache/ en zooi
