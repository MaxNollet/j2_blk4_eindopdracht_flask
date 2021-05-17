# insert_genepanel(tuple(str)): bool
import os
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from gaps.models import Gene
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

    def __init__(self, engine, session):
        self.engine = engine
        self.session = session
        self.test()

    def create__engine(self):
        self.engine = create_engine(
            "postgresql://maxn:blaat1234@bio-inf.han.nl:5432/maxn", echo=True)

    def create__session(self):
        print("Creating session")
        Session = sessionmaker(bind=self.engine)
        self.session = Session()

    def test(self):
        gene = Gene(id=None, ncbi_gene_id='8139', hgnc_symbol='GAAAN',
                    in_genepanel=True)
        self.session.add(gene)
        self.session.commit()


def main():
    # insertGene()
    # createEngine()
    # engine = create_engine(
    #     "postgresql://maxn:blaat1234@bio-inf.han.nl:5432/maxn", echo=True)
    # print("Creating session")
    # Session = sessionmaker(bind=engine)
    # session = Session()
    # Alchemy(engine, session)
    Alchemy(engine, session)


main()
