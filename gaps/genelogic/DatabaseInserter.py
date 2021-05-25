# insert_genepanel(tuple(str)): bool
import bdb
import os
from gaps.genelogic import reader
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from gaps.models import *
from sqlalchemy.exc import IntegrityError

import psycopg2


# https://docs.sqlalchemy.org/en/14/orm/tutorial.html


# def insertGene():
#     pass
#
# def session():
#     # Session = sessionmaker(bind=engine)
#     # Session =
#     pass
#
# def createEngine():
#     postgre = os.environ.get('SQLALCHEMY_DATABASE_URI')
#     postgre = os.getenv("SQLALCHEMY_DATABASE_URI")
#
#     print(postgre)
#
#
#     Session = sessionmaker(bind=engine)
#     print("Creating session")
#     session = Session()
#     t = Gene(id=None, ncbi_gene_id='8139', hgnc_symbol='GAAN', in_genepanel=True)
#     # session.add(t)
#
#     tt = session.query(Gene).filter_by(ncbi_gene_id='8139').count()
#     session.commit()
#     print(tt)
#
#
#     session.close()

# for row in session.query(Gene, Gene.ncbi_gene_id).all():
#     print(row.Gene, row.ncbi_gene_id)
class Alchemy:

    def __init__(self):
        self.engine = None
        self.session = None

    def create__engine(self):
        self.engine = create_engine()

    def __create__engine(self, db):
        self.engine = create_engine(db, echo=True)

    def __create__session(self):
        print("Creating session")
        Session = sessionmaker(bind=self.engine)
        self.session = Session()

    def test(self):
        self.__create__engine(
            "postgresql://maxn:blaat1234@bio-inf.han.nl:5432/maxn")
        self.__create__session()
        gene = Gene(id=None, ncbi_gene_id='8139', hgnc_symbol='GAPS',
                    in_genepanel=True)
        self.session.add(gene)
        self.session.commit()


def update_genepanel_v2():
    """A function which inserts every line of a file
       into the database.
    """
    print("Triggered")
    path = os.path.join(os.path.dirname(__file__), "GenPanelOverzicht_DG-3.1.0_HAN_original_tsv.txt")
    if os.path.exists(path):
        file = reader.get_reader(path)
        session = db.session
        try:
            cached_genes = dict((gene.hgnc_symbol, gene) for gene in Gene.query.all())
            cached_aliases = dict((alias.hgnc_symbol, alias) for alias in Alias.query.all())
            cached_genepanel_symbols = dict((genepanel_symbol.symbol, genepanel_symbol) for genepanel_symbol in GenepanelSymbol.query.all())
            for line in file:
                gene = line.gene
                if gene.hgnc_symbol in cached_genes:
                    gene = cached_genes[gene.hgnc_symbol]
                else:
                    cached_genes[gene.hgnc_symbol] = gene
                # Insert aliases.
                for alias in line.alias:
                    if alias.hgnc_symbol in cached_aliases:
                        gene.aliases.append(cached_aliases[alias.hgnc_symbol])
                    else:
                        gene.aliases.append(alias)
                        cached_aliases[alias.hgnc_symbol] = alias
                # # Insert genepanel symbols.
                # genepanel_symbol = GenepanelSymbol(symbol=line.p_symbol.symbol)
                # if genepanel_symbol.symbol in cached_genepanel_symbols:
                #     gene.genepanel_symbol = cached_genepanel_symbols[genepanel_symbol.symbol]
                # else:
                #     gene.genepanel_symbol = genepanel_symbol
                #     cached_genepanel_symbols[genepanel_symbol.symbol] = genepanel_symbol

                session.add(gene)
                session.flush()
            session.commit()
        except IntegrityError:
            session.rollback()






        # count = len(file)
        # counter = 1
        # for line in file:
        #     print(f"\r{counter} / {count}")
        #     insert = line.gene
        #     for alias in line.alias:
        #         retrieved = Alias.query.filter(Alias.hgnc_symbol == alias.hgnc_symbol).first()
        #         if retrieved is None:
        #             insert.aliases.append(alias)
        #         else:
        #             insert.aliases.append(retrieved)
        #     db.session.add(insert)
        #     db.session.flush()
        #     counter += 1
        # db.session.commit()
            # for alias in line.alias:
            #     print(db.session.get(Alias, alias.hgnc_symbol))
            # # [insert.aliases.append(alias) for alias in line.alias]
            # # db.session.add(insert)
            # db.session.merge(insert)
            # db.session.flush()
        # db.session.commit()


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

                            for test in testen:  # test = [InheritanceType(id=None, type='AD')]
                                print(k, "Genepanel")
                                duplicate = test.query.filter_by(
                                    type=test).first()
                                print(duplicate, " duplicate")
                                if duplicate is None:
                                    k.inheritance_types.append(test)
                                    db.session.add(k)
                                    db.session.commit()
                                else:
                                    k.inheritance_types.append(duplicate.id)
                                # print(k, "Genepanel")
                                # duplicate = k.query.filter_by(
                                #     abbreviation=k.abbreviation).first()
                                # print(duplicate, " duplicate")
                                # if duplicate is None:
                                #     dp = test.query.filter_by(
                                #         type=test.type).first()
                                #     if dp is None:
                                #         k.inheritance_types.append(test)
                                #         db.session.add(k)
                                #         db.session.commit()
                                #     else:
                                #         k.inheritance_types.append(test.id)
                                #         db.session.add(k)
                                #         db.session.commit()
                                # else:
                                #     print("else check")
                                #     k.inheritance_types.append(duplicate.id)
                                #     db.session.add(k)
                                #     db.session.commit()
                            testen = []
                        if isinstance(k,
                                      InheritanceType):  # kijkt of het het juiste object is
                            # ih = k
                            testen.append(k)
                            # >> > peter = User.query.filter_by(
                            #     username='peter').first()
                            # moet unique zijn
                            # duplicate = InheritanceType.query.filter_by(
                            #     type=k.type).first()
                            # if duplicate is None:  # zodat er geen dubbele in komen
                            #     print("=-=-=", k)
                            #     line.gene.genepanels.append(k)
                            #     # de tussen tabelen moetne nog
                            #     db.session.add(k)
                            #     db.session.commit()
            db.session.commit()
        print("Done")
    else:
        print("No file!")

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
