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


def updateGenpanel():
    # Pad opvragen waar dit bestand staat waar we nu in zitten,
    # vervolgens naam van tsv-bestand eraan toevoegen.
    path = os.path.join(os.path.dirname(__file__), "GenPanelOverzicht_DG-3.1.0_HAN_original_tsv.txt")
    # path = "/Users/lean/Documenten/School/Flask/Course8_project/gaps/genelogic/GenPanelOverzicht_DG-3.1.0_HAN_original_tsv.txt"
    if os.path.exists(path):
        data = reader.get_reader(path)
        # print(len(data))

        # gene = Gene(id=None, ncbi_gene_id='8139', hgnc_symbol='GAAPS',
        #             in_genepanel=True)
        # db importeer je al met 'from gaps.models import *'.
        # Note: het werkte eerst niet omdat je deze functie deed aanroepen in
        # __init__.py terwijl de app nog niet volledig was opgestart, ook
        # was de database nog niet geïnitialiseerd toen je de aanroep deed.
        # Oplossing: deze methode in /query van blueprint_query aanroepen
        # nadat de app volledig is opgestart en de database is geïnitialiseerd.
        # db.session.add(gene)
        # db.session.delete(gene)  # hij loopt hier vast
        # db.session.commit()  # https://flask-sqlalchemy.palletsprojects.com/en/2.x/queries/

        for line in data:
            # print("iets")
            # g = line.gene
            # gene_tabel = Gene(ncbi_gene_id=line.gene.ncbi_gene_id,
            #                   hgnc_symbol=line.gene.hgnc_symbol,
            #                   in_genepanel=line.gene.in_genepanel)
            # print(gene_tabel)
            # # db.session.add(g)
            # db.session.add(gene_tabel)
            # db.session.commit()
            # # gene_id = g.id
            # current_gene = Gene.query.filter_by(ncbi_gene_id=line.gene.ncbi_gene_id).first()
            # print(current_gene.id, " werkt dit?")
            # # db.session.commit()
            # p = line.p_symbol       # Genepanel Symbol
            # p = GenepanelSymbol(symbol=p.symbol, gene_id=current_gene.id)
            # # p.gene = gene_id
            # # print(gene_id, p)
            # db.session.add(p)
            # db.session.commit()
            if len(line.alias) >= 1:
                for i in line.alias:  # loop over de meerdere aliases
                    print(i, " ja")
                    # t = Alias()
                    # t.genes.append(i)
                    # kaas = Alias(hgnc_symbol=i.hgnc_symbol)
                    # db.session.add(kaas)
                    # db.session.commit() # je kan niet het id ophalen
                    # want hij zit nog niet in de db (None dus), maar
                    # als je het Alias er inzit waardoor er een id bij
                    # komt is het Key (hgnc_symbol)=(A2MD) already exists.?!
                    # i.genes.append(line.gene.id)
                    # db.session.commit()
                    # p = i.genes.append(line.gene.id)  probeersel
                    # i.genes.append(line.gene.id) # sqlalchemy.exc.IntegrityError: (psycopg2.errors.UniqueViolation) duplicate key value violates unique constraint "alias_symbol_unique"
                    # DETAIL:  Key (hgnc_symbol)=(A2MD) already exists.
                    #                     i.genes.append(current_gene.id)     ((c is not None) and instance_state(c) or None, c)
                    # AttributeError: 'int' object has no attribute '_sa_instance_state'

                    i.genes.append(line.gene.id)
                    # https://stackoverflow.com/questions/25668092/flask-sqlalchemy-many-to-many-insert-data
                    db.session.add(i)
                    db.session.commit()
                    # Alias.id.append(gene_id)
                    # genes moet nu gevuld worden
                    # db.session.add(t)
                    # db.session.commit()
            if len(line.panel) >= 1:
                # https://www.tutorialspoint.com/sqlalchemy/sqlalchemy_orm_many_to_many_relationships.htm
                # https://flask-sqlalchemy.palletsprojects.com/en/2.x/queries/
                for j in line.panel:
                    testen = []
                    for k in j:
                        print(testen, k, " testen, k")
                        # print(isinstance(k,
                        #                  InheritanceType))  # K = Inheritance en GenePanel
                        # if k.
                        if isinstance(k, Genepanel):
                            pass

# def main():
#     path = "/Users/lean/Documenten/School/Flask/Course8_project/gaps/genelogic/GenPanelOverzicht_DG-3.1.0_HAN_original_tsv.txt"
# #
#     data = reader.get_reader(path)
#     print(data[0])
#     for line in data:
#         print(line.gene)
#         t = line.gene
#
#         print(line.p_symbol)
#         if len(line.alias) >= 1:
#             for i in line.alias:
#                 print("***", i)
#         if len(line.panel) >= 1:
#             for j in line.panel:
#                 for k in j:
#                     print("=-=", k)


#     # insertGene()
#     # createEngine()
#     # engine = create_engine(
#     #     "postgresql://maxn:blaat1234@bio-inf.han.nl:5432/maxn", echo=True)
#     # print("Creating session")
#     # Session = sessionmaker(bind=engine)
#     # session = Session()
#     # Alchemy(engine, session)
#     # t = Alchemy()
#     # t.updateGenpanel()
#     gene = Gene(id=None, ncbi_gene_id='8139', hgnc_symbol='GAPS',
#                     in_genepanel=True)
#     db.session.add(gene)
#     SQLALCHEMY_DATABASE_URI = os.environ.get("SQLALCHEMY_DATABASE_URI")
#     print(SQLALCHEMY_DATABASE_URI)


# main()
