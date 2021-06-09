from flask import Blueprint, render_template

blueprint_help = Blueprint("blueprint_help", __name__)


@blueprint_help.route("/help")
def route_help():
    """A function which handles requests to the '/help'-
       route for the webapp.

    :return Rendered template of 'template_help.html'.
    """
    return render_template("template_help.html", active_nav="help")
