from sqlalchemy import update, bindparam

from gaps.genelogic.statement_groups import statement_group
from gaps.models import Gene


@staticmethod
@statement_group(table="genepanel_symbol_id")
def _relation_genepanel_symbol_id():
    return update(Gene).where(
        Gene.hgnc_symbol == bindparam("hgnc_symbol")
    ).values(
        {Gene.genepanel_symbol_id: bindparam("genepanel_symbol_id")}
    )
