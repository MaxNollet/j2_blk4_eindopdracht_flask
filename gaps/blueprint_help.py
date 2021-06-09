from flask import Blueprint, render_template

blueprint_help = Blueprint("blueprint_help", __name__)


@blueprint_help.route("/help")
def route_help():
    return render_template("template_help.html", active_nav="help")
