# insert_genepanel(tuple(str)): bool
import os
import time
from typing import List, Mapping

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

    def __init__(self):
        self.session: Session = db.session

    def insert_genepanel(self, rows: List[Mapping[str, Mapping[str, str]]]):
        genes: List[dict] = list()
        aliases: List[dict] = list()
        genepanel_symbols: List[dict] = list()
        for row in rows:
            gene = row.get("gene")
            if gene is not None and len(gene) > 0:
                genes.append(gene)
            alias = row.get("aliases")
            if alias is not None and len(alias) > 0:
                aliases.extend(alias)
            genepanel_symbol = row.get("genepanel_symbol")
            if genepanel_symbol is not None:
                genepanel_symbols.append(genepanel_symbol)
        starttijd = time.perf_counter()

        gene_ids = self.insert(insert_stmt=self.stmt_gene_insert, values=genes, return_ids=True,
                               select_stmt=self.stmt_gene_select, col_as_key="hgnc_symbol")
        alias_ids = self.insert(insert_stmt=self.stmt_alias_insert, values=aliases, return_ids=True,
                                select_stmt=self.stmt_alias_select, col_as_key="hgnc_symbol")
        genepanel_symbols_ids = self.insert(insert_stmt=self.stmt_genepanel_symbol_insert, values=genepanel_symbols,
                                            return_ids=True, select_stmt=self.stmt_genepanel_symbol_select,
                                            col_as_key="symbol")
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
            # test = {result[0]: result[1] for result in results}
            return {result[0]: result[1] for result in results}

    @staticmethod
    def link_gene_alias(session: Session, rows: List[dict]):
        gene_symbols = list()
        alias_symbols = list()
        for row in rows:
            gene_symbols.append(row.get("gene").get("hgnc_symbol"))
            for alias in row.get("aliases"):
                alias_symbols.append(alias.get("hgnc_symbol"))

        gene_objects = dict(
            (gene.hgnc_symbol, gene) for gene in Gene.query.filter(Gene.hgnc_symbol.in_(gene_symbols)).all())
        alias_objects = dict(
            (alias.hgnc_symbol, alias) for alias in Alias.query.filter(Alias.hgnc_symbol.in_(alias_symbols)).all())

        combos = list()
        for row in rows:
            id_gene = gene_objects.get(row.get("gene").get("hgnc_symbol")).id
            for alias in row.get("aliases"):
                id_alias = alias_objects.get(alias.get("hgnc_symbol")).id
                combos.append({"gene_id": id_gene, "alias_id": id_alias})
        statement = insert(t_gene_alias).on_conflict_do_nothing()
        session.execute(statement=statement, params=combos)
        # session.commit()


def update_genepanel_v2():
    """A function which inserts every line of a file
       into the database.
    """
    print("Triggered")
    path = os.path.join(os.path.dirname(__file__), "GenPanelOverzicht_DG-3.1.0_HAN_original_tsv.txt")
    if os.path.exists(path):
        file = reader.get_reader(path)
        session = db.session
        all_rows: List[Mapping[str, Mapping[str, str]]] = list()
        for line in file:
            gene = line.gene
            row = {"gene": {"ncbi_gene_id": gene.ncbi_gene_id, "hgnc_symbol": gene.hgnc_symbol, "in_genepanel": True},
                   "aliases": []}
            for alias in line.alias:
                if alias.hgnc_symbol is not None and alias.hgnc_symbol != "":
                    row["aliases"].append({"hgnc_symbol": alias.hgnc_symbol})
            row["genepanel_symbol"] = {"symbol": line.p_symbol.symbol}
            all_rows.append(row)

        print("Executing . . .")
        inserter = DatabaseInserter()
        inserter.insert_genepanel(all_rows)
        print("Done")
