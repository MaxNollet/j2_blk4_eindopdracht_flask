# insert_genepanel(tuple(str)): bool
import os
import time
from typing import List, Mapping, Tuple

from sqlalchemy import select, bindparam
from sqlalchemy.dialects.postgresql import insert
from sqlalchemy.orm import Session

from gaps.genelogic import reader
from gaps.models import *


class DatabaseInserter:
    stmt_gene_insert = insert(Gene).on_conflict_do_nothing()
    stmt_gene_select = select(Gene.hgnc_symbol, Gene.id).where(Gene.hgnc_symbol.in_(bindparam("values")))
    stmt_alias_insert = insert(Alias).on_conflict_do_nothing()
    stmt_alias_select = select(Alias.hgnc_symbol, Alias.id).where(Alias.hgnc_symbol.in_(bindparam("values")))
    stmt_genepanel_symbol_insert = insert(GenepanelSymbol).on_conflict_do_nothing()
    stmt_genepanel_symbol_select = select(GenepanelSymbol.symbol, GenepanelSymbol.id).where(
        GenepanelSymbol.symbol.in_(bindparam("values")))
    stmt_genepanel_insert = insert(Genepanel).on_conflict_do_nothing()
    stmt_genepanel_select = select(Genepanel.abbreviation, Genepanel.id).where(
        Genepanel.abbreviation.in_(bindparam("values")))
    stmt_inheritance_type_insert = insert(InheritanceType).on_conflict_do_nothing()
    stmt_inheritance_type_select = select(InheritanceType.type, InheritanceType.id).where(
        InheritanceType.type.in_(bindparam("values")))

    stmt_relation_gene_alias = insert(t_gene_alias).on_conflict_do_nothing()

    def __init__(self):
        self.session: Session = db.session

    def insert_genepanel(self, file):
        # New improved data structures for the db.
        all_genes = list()
        all_aliases = list()
        all_genepanel_symbols = list()
        all_genepanels = list()
        all_inheritace_types = list()
        relation_gene_alias = list()
        relation_gene_symbol = list()
        relation_gene_genepanel = list()
        relation_genepanel_inheritance = list()
        # Extract variables and seperate into categories.
        for line in file:
            all_genes.append(
                {"ncbi_gene_id": line.gene.ncbi_gene_id, "hgnc_symbol": line.gene.hgnc_symbol, "in_genepanel": True})
            for alias in line.alias:
                if alias.hgnc_symbol is not None and alias.hgnc_symbol != "":
                    all_aliases.append({"hgnc_symbol": alias.hgnc_symbol})
                    relation_gene_alias.append({"gene_id": line.gene.hgnc_symbol, "alias_id": alias.hgnc_symbol})
            all_genepanel_symbols.append({"symbol": line.p_symbol.symbol})
            for panel in line.panel:
                [all_inheritace_types.append({"type": element.type}) for element in panel[:-1]]
                all_genepanels.append({"abbreviation": panel[-1].abbreviation})
        # Generate combinations

        starttijd = time.perf_counter()

        gene_ids = self.insert(insert_stmt=self.stmt_gene_insert, values=all_genes, return_ids=True,
                               select_stmt=self.stmt_gene_select, col_as_key="hgnc_symbol")
        alias_ids = self.insert(insert_stmt=self.stmt_alias_insert, values=all_aliases, return_ids=True,
                                select_stmt=self.stmt_alias_select, col_as_key="hgnc_symbol")
        genepanel_symbols_ids = self.insert(insert_stmt=self.stmt_genepanel_symbol_insert, values=all_genepanel_symbols,
                                            return_ids=True, select_stmt=self.stmt_genepanel_symbol_select,
                                            col_as_key="symbol")
        genepanel_ids = self.insert(insert_stmt=self.stmt_genepanel_insert, values=all_genepanels, return_ids=True,
                                    select_stmt=self.stmt_genepanel_select, col_as_key="abbreviation")
        inheritance_ids = self.insert(insert_stmt=self.stmt_inheritance_type_insert, values=all_inheritace_types,
                                      return_ids=True, select_stmt=self.stmt_inheritance_type_select,
                                      col_as_key="type")
        # Opgehaalde primary keys gebruiken om relatie te updaten.
        test = self.generate_combinations(relation_gene_alias, {"gene_id": gene_ids, "alias_id": alias_ids})
        self.relationship(self.stmt_relation_gene_alias, test)
        self.session.commit()
        # self.session.commit()
        # DatabaseInserter.insert_genepanel_symbol(session, genepanel_symbols)
        # DatabaseInserter.link_gene_alias(session, rows)
        print(f"Verwerktijd: {time.perf_counter() - starttijd}")

    def insert(self, insert_stmt: insert, values: List[dict], return_ids: bool = False,
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

    def relationship(self, insert_stmt: insert, values: Tuple[dict]):
        self.session.execute(statement=insert_stmt, params=values)

    def generate_combinations(self, original_combinations: List[dict],
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
