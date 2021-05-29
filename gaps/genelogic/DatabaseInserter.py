# insert_genepanel(tuple(str)): bool
import os
from typing import List, Mapping

import time

from sqlalchemy import create_engine, bindparam
from sqlalchemy.dialects.postgresql import insert
from sqlalchemy import select
from sqlalchemy.exc import IntegrityError
from sqlalchemy.orm import sessionmaker, Session

from gaps.genelogic import reader
from gaps.models import *


class DatabaseInserter:
    @staticmethod
    def insert_genepanel(session: Session, rows: List[Mapping[str, Mapping[str, str]]]):
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
        DatabaseInserter.insert_genes(session, genes)
        DatabaseInserter.insert_aliases(session, aliases)
        DatabaseInserter.insert_genepanel_symbol(session, genepanel_symbols)
        DatabaseInserter.link_gene_alias(session, rows)

        print(f"Verwerktijd: {time.perf_counter() - starttijd}")

    @staticmethod
    def insert_genes(session: Session, genes: List[dict]):
        """Insert genes into the database."""
        statement = insert(Gene).on_conflict_do_nothing()
        session.execute(statement=statement, params=genes)
        session.commit()

    @staticmethod
    def insert_aliases(session: Session, aliases: List[dict]):
        """Insert aliases into the database."""
        statement = insert(Alias).on_conflict_do_nothing()
        session.execute(statement=statement, params=aliases)
        session.commit()

    @staticmethod
    def insert_genepanel_symbol(session: Session, symbols: List[dict]):
        """Insert genepanel symbols into the database."""
        statement = insert(GenepanelSymbol).on_conflict_do_nothing()
        session.execute(statement=statement, params=symbols)
        session.commit()


    @staticmethod
    def link_gene_alias(session: Session, rows: List[dict]):
        gene_symbols = list()
        alias_symbols = list()
        for row in rows:
            gene_symbols.append(row.get("gene").get("hgnc_symbol"))
            for alias in row.get("aliases"):
                alias_symbols.append(alias.get("hgnc_symbol"))

        gene_objects = dict((gene.hgnc_symbol, gene) for gene in Gene.query.filter(Gene.hgnc_symbol.in_(gene_symbols)).all())
        alias_objects = dict((alias.hgnc_symbol, alias) for alias in Alias.query.filter(Alias.hgnc_symbol.in_(alias_symbols)).all())

        combos = list()
        for row in rows:
            id_gene = gene_objects.get(row.get("gene").get("hgnc_symbol")).id
            for alias in row.get("aliases"):
                id_alias = alias_objects.get(alias.get("hgnc_symbol")).id
                combos.append({"gene_id": id_gene, "alias_id": id_alias})
        statement = insert(t_gene_alias).on_conflict_do_nothing()
        session.execute(statement=statement, params=combos)
        session.commit()


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
        DatabaseInserter.insert_genepanel(session, all_rows)
        print("Done")
