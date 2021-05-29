import os
import time
from typing import List, Mapping, Tuple

from sqlalchemy import select
from sqlalchemy.dialects.postgresql import insert
from sqlalchemy.orm import Session

from gaps.genelogic import reader, SelectStatements, InsertStatements
from gaps.models import db


class DatabaseInserter(SelectStatements, InsertStatements):
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

        gene_ids = self.insert_values(self._insert_gene(), all_genes, True, self._select_gene(), "hgnc_symbol")
        alias_ids = self.insert_values(self._insert_alias(), all_aliases, True, self._select_alias(), "hgnc_symbol")
        genepanel_symbols_ids = self.insert_values(self._insert_genepanel_symbol(), all_genepanel_symbols, True,
                                                   self._select_genepanel_symbol(), "symbol")
        genepanel_ids = self.insert_values(self._insert_genepanel(), all_genepanels, True, self._select_genepanel(),
                                    "abbreviation")
        inheritance_ids = self.insert_values(self._insert_inheritance_type(), all_inheritace_types, True,
                                             self._select_inheritance_type(), "type")
        # Opgehaalde primary keys gebruiken om relatie te updaten.
        pks_gene_alias = self.combine(relation_gene_alias, {"gene_id": gene_ids, "alias_id": alias_ids})
        pks_gene_genepanel = self.combine(relation_gene_genepanel, {"gene_id": gene_ids, "genepanel_id": genepanel_ids})
        pks_genepanel_inheritance = self.combine(relation_genepanel_inheritance, {"genepanel_id": genepanel_ids,
                                                                                  "inheritance_type_id": inheritance_ids})
        self.insert_relationship(self._insert_relation_gene_alias(), pks_gene_alias)
        self.insert_relationship(self._insert_relation_gene_genepanel(), pks_gene_genepanel)
        self.insert_relationship(self._insert_relation_genepanel_inheritance(), pks_genepanel_inheritance)
        self.session.commit()
        print(f"Verwerktijd: {time.perf_counter() - starttijd}")

    def insert_values(self, insert_stmt: insert, values: List[dict], return_ids: bool = False,
                      select_stmt: select = None, col_as_key: str = None) -> Mapping[str, int]:
        """A method which inserts values into the database and
           returns the primary keys of all inserted values if
           return_ids is set to True.

        :param insert_stmt Statement to use when inserting values (Insert).
        :param values All values to be inserted into the database (List[dict]).
        :param return_ids Indication to return inserted ID's or not (Bool).
        :param select_stmt Statement to use when selecting values (Select).
        :param col_as_key Column to use for the keys in the returned dictionary (str).
        :return All inserted primary keys bound to corresponding values (Mapping[str, int]).
        """
        self.session.execute(statement=insert_stmt, params=values)
        if return_ids is True and select_stmt is not None and col_as_key is not None:
            values = {"values": [value[col_as_key] for value in values]}
            results = self.session.execute(statement=select_stmt, params=values)
            return {result[0]: result[1] for result in results}

    def insert_relationship(self, insert_stmt: insert_values, values: Tuple[dict]):
        self.session.execute(statement=insert_stmt, params=values)

    @staticmethod
    def combine(original_combinations: List[dict],
                primary_keys: Mapping[str, Mapping[str, int]]) -> Tuple[dict]:
        """A method which updates combinations between two tables
           using primary keys retrieved from the database.

        :param original_combinations Combinations between two tables using normal values (List[dict]).
        :param primary_keys Primary keys to use when updating normal values (List[dict]).
        :return combinations of primary keys to insert into in-between-tables (List[Dict]).
        """
        keys = primary_keys.keys()
        for combo in original_combinations:
            for key in keys:
                original_value = combo[key]
                combo[key] = primary_keys[key][original_value]
        return tuple(original_combinations)


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
