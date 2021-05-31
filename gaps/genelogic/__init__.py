from gaps.gene_searcher import GeneSearcher_code as gs
from gaps.genelogic.statement_groups import statement_group
from gaps.genelogic.select_statements import SelectStatements
from gaps.genelogic.insert_statements import InsertStatements


def search_query(query):
    results = gs.results_query(query)

    return results
