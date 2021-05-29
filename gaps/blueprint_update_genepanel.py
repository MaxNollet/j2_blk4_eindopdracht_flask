from flask import Blueprint, render_template

from gaps.models import db, Gene, Alias, Genepanel, t_gene_alias
from sqlalchemy import func

from gaps.genelogic.DatabaseInserter import update_genepanel_v2

blueprint_update_genepanel = Blueprint("blueprint_update_genepanel", __name__)


@blueprint_update_genepanel.route("/update_genepanel")
def update_genepanel():
    """A function which handles requests for the '/update_genepanel'
       -route for the webapp.
    """
    update_genepanel_v2()
    # test = Alias.query.join(t_gene_alias).join(Alias).filter(t_gene_alias.c.alias_id == Alias.id).all()
    # print(test)
    # print(db.session.query(func.count(Gene.id)).group_by(Gene.genepanels).all())
    # most_occuring_alias = db.session.query(func.count(Alias.id)).group_by(Alias.genes).all()
    # print(f"Genes: {count_genes}")
    # print(f"Aliases: {count_aliases}")
    # print(f"Genepanels: {count_genepanels}")
    return render_template("template_update_genepanel.html",
                           active_nav="genepanel",
                           unique_genes=Gene.query.count(),
                           unique_aliases=Alias.query.count(),
                           unique_genepanels=Genepanel.query.count())
