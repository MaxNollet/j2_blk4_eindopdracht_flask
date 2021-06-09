from flask import Blueprint, render_template
from sqlalchemy.exc import OperationalError

from gaps.models import Gene, Alias, Genepanel

blueprint_update_genepanel = Blueprint("blueprint_update_genepanel", __name__)


@blueprint_update_genepanel.route("/update_genepanel")
def update_genepanel():
    """A function which handles requests for the '/update_genepanel'-
       route for the webapp.

    :return Rendered template of 'template_update_genepanel.html'.
    """
    try:
        genes = Gene.query.count()
        aliases = Alias.query.count()
        genepanels = Genepanel.query.count()
        return render_template("template_update_genepanel.html",
                               active_nav="genepanel",
                               unique_genes=genes,
                               unique_aliases=aliases,
                               unique_genepanels=genepanels)
    except OperationalError:
        return render_template("template_update_genepanel.html",
                               active_nav="genepanel")
