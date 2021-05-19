# insert_genepanel(tuple(str)): bool
import os
from gaps.genelogic import reader
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from gaps.models import *

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
#     engine = create_engine("postgresql://maxn:blaat1234@bio-inf.han.nl:5432/maxn", echo=True)
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


def updateGenpanel(db):
    print("test")
    path = "/Users/lean/Documenten/School/Flask/Course8_project/gaps/genelogic/GenPanelOverzicht_DG-3.1.0_HAN_original_tsv.txt"
    if os.path.exists(path):
        data = reader.get_reader(path)
        # print(len(data))
        # gene = Gene(id=None, ncbi_gene_id='8139', hgnc_symbol='GAAPS',
        #             in_genepanel=True)
        # db.session.add(gene)
        # db.session.commit()
        for line in data:
            print("iets")
            g = line.gene
            db.session.add(g)
            db.session.commit()
            p = line.p_symbol
            db.session.add(p)
            if len(line.alias) >= 1:
                for i in line.alias:
                    db.session.add(i)
            if len(line.panel) >= 1:
                for j in line.panel:
                    for k in j:
                        db.session.add(k)
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
