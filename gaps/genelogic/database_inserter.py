import os
import time
from typing import List, Mapping, Tuple

from sqlalchemy.orm import Session

from gaps.genelogic import reader, SelectStatements, InsertStatements, StatementGroups
from gaps.genelogic.statement_groups import InsertStatementNotDefined, SelectStatementNotDefined,\
    ColumnAsKeyNotDefined, StatementGroupNotDefined
from gaps.models import db


class DatabaseInserter(StatementGroups, SelectStatements, InsertStatements):
    """A class which contains several methods for inserting
       values into the database. This class extends the classes
       SelectStatements and InsertStatements to rely on statements
       defined in thode classes.
    """

    def __init__(self, session: Session = None):
        """Constructor of the object. Uses the session-object
           provided to communicate with the database or creates
           its own session if no session is provided.

        :param session Session-object used to communicate with the database (Session).
        """
        if session is None:
            self.session: Session = db.session
        else:
            self.session = session
        self.ids = dict()

    def insert_genepanel(self, file):
        # New improved data structures for the db.
        all_genes = list()
        all_aliases = list()
        all_genepanel_symbols = list()
        all_genepanels = list()
        all_inheritace_types = list()
        relation_gene_alias = list()
        relation_gene_genepanel = list()
        relation_genepanel_inheritance = list()
        # Extract variables and separate into categories.
        for line in file:
            all_genes.append(
                {"ncbi_gene_id": line.gene.ncbi_gene_id, "hgnc_symbol": line.gene.hgnc_symbol, "in_genepanel": True})
            for alias in line.alias:
                if alias.hgnc_symbol is not None and alias.hgnc_symbol != "":
                    all_aliases.append({"hgnc_symbol": alias.hgnc_symbol})
                    relation_gene_alias.append({"gene_id": line.gene.hgnc_symbol, "alias_id": alias.hgnc_symbol})
            all_genepanel_symbols.append({"symbol": line.p_symbol.symbol})
            for panel in line.panel:
                all_genepanels.append({"abbreviation": panel[-1].abbreviation})
                [all_inheritace_types.append({"type": element.type}) for element in panel[:-1]]
                relation_gene_genepanel.append(
                    {"gene_id": line.gene.hgnc_symbol, "genepanel_id": panel[-1].abbreviation})
                [relation_genepanel_inheritance.append(
                    {"genepanel_id": panel[-1].abbreviation, "inheritance_type_id": element.type}) for element in
                    panel[:-1]]
        # Generate combinations

        starttijd = time.perf_counter()

        self.ids["gene_id"] = self.insert_values("gene", all_genes, True)
        self.ids["alias_id"] = self.insert_values("alias", all_aliases, True)
        self.ids["genepanel_symbol_id"] = self.insert_values("genepanel_symbol", all_genepanel_symbols, True)
        self.ids["genepanel_id"] = self.insert_values("genepanel", all_genepanels, True)
        self.ids["inheritance_type_id"] = self.insert_values("inheritance_type", all_inheritace_types, True)
        # Opgehaalde primary keys gebruiken om relatie te updaten.
        pks_gene_alias = self.combine(relation_gene_alias, ("gene_id", "alias_id"))
        pks_gene_genepanel = self.combine(relation_gene_genepanel, ("gene_id", "genepanel_id"))
        pks_genepanel_inheritance = self.combine(relation_genepanel_inheritance,
                                                 ("genepanel_id", "inheritance_type_id"))
        self.insert_values("gene_alias", pks_gene_alias)
        self.insert_values("genepanel_gene", pks_gene_genepanel)
        self.insert_values("genepanel_inheritance", pks_genepanel_inheritance)
        # self.session.commit()
        print(f"Verwerktijd: {time.perf_counter() - starttijd}")

    def insert_values(self, table_name: str, values: List[dict], return_ids: bool = False) -> Mapping[str, int]:
        """A method which inserts values into the database and
           returns the primary keys of all inserted values if
           return_ids is set to True.

        :param table_name Name of the table to insert values into (str).
        :param values All values to be inserted into the database (List[dict]).
        :param return_ids Indication to return inserted ID's or not (Bool).
        :return All inserted primary keys bound to corresponding values (Mapping[str, int]).
        """
        stmt_group = self._statement_groups.get(table_name)
        if stmt_group:
            insert_stmt = stmt_group.get("insert")
            if insert_stmt:
                self.session.execute(statement=insert_stmt, params=values)
            else:
                raise InsertStatementNotDefined(table_name)
            if return_ids is True:
                col_as_key = stmt_group.get("column_as_key")
                select_stmt = stmt_group.get("select")
                if col_as_key and select_stmt:
                    values = {"values": [value[col_as_key] for value in values]}
                    results = self.session.execute(statement=select_stmt, params=values)
                    return {result[0]: result[1] for result in results}
                else:
                    if not select_stmt:
                        raise SelectStatementNotDefined(table_name)
                    if not col_as_key:
                        raise ColumnAsKeyNotDefined(table_name)
        else:
            raise StatementGroupNotDefined(table_name)

    def combine(self, original_combinations: List[dict],
                primary_keys: Tuple[str, str]) -> List[dict]:
        """A method which updates combinations between two tables
           using primary keys retrieved from the database.

        :param original_combinations Combinations between two tables using normal values (List[dict]).
        :param primary_keys Primary keys to use when updating normal values (List[dict]).
        :return combinations of primary keys to insert into in-between-tables (List[Dict]).
        """
        for combo in original_combinations:
            for key in primary_keys:
                original_value = combo[key]
                combo[key] = self.ids[key][original_value]
        return original_combinations


def update_genepanel_v2():
    """A function which inserts every line of a file
       into the database.
    """
    print("Triggered")
    path = os.path.join(os.path.dirname(__file__), "GenPanelOverzicht_DG-3.1.0_HAN_original_tsv.txt")
    if os.path.exists(path):
        file = reader.get_reader(path)
        print("Executing . . .")
        inserter = DatabaseInserter()
        inserter.insert_genepanel(file)
        print("Done")
