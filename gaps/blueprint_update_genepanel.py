from flask import Blueprint, render_template

blueprint_update_genepanel = Blueprint("blueprint_update_genepanel", __name__)


@blueprint_update_genepanel.route("/update_genepanel")
def update_genepanel():
    """A function which handles requests for the '/update_genepanel'
       -route for the webapp.
    """
    return render_template("template_update_genepanel.html", active_nav="genepanel")
