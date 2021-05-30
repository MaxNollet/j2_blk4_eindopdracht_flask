from gaps.gene_searcher import GeneSearcher_code as gs
from gaps.genelogic.insert_statements import InsertStatements
from gaps.genelogic.select_statements import SelectStatements
from gaps.genelogic.statement_groups import StatementGroups


def search_query(query):
    results = gs.results_query(query)

    return results
