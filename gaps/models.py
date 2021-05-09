# coding: utf-8
from sqlalchemy import Boolean, Column, Date, ForeignKey, Integer, String, Table, text
from sqlalchemy.orm import relationship
from flask_sqlalchemy import SQLAlchemy
from sqlalchemy.ext.declarative import declarative_base

db = SQLAlchemy()

Base = declarative_base()
metadata = Base.metadata


class Gene(Base):
    __tablename__ = 'gene'
    __table_args__ = {'schema': 'eindopdracht'}

    id = Column(Integer, primary_key=True, server_default=text("nextval('eindopdracht.gene_id_seq'::regclass)"))
    ncbi_gene_id = Column(String(60), nullable=False)
    hgnc_symbol = Column(String(15), nullable=False)
    in_genepanel = Column(Boolean, nullable=False)

    genepanels = relationship('Genepanel', secondary='eindopdracht.genepanel_gene')


class Genepanel(Base):
    __tablename__ = 'genepanel'
    __table_args__ = {'schema': 'eindopdracht'}

    id = Column(Integer, primary_key=True, server_default=text("nextval('eindopdracht.genepanel_id_seq'::regclass)"))
    afkorting = Column(String(40), nullable=False)
    naam = Column(String(100), nullable=False)


class Journal(Base):
    __tablename__ = 'journal'
    __table_args__ = {'schema': 'eindopdracht'}

    id = Column(Integer, primary_key=True, server_default=text("nextval('eindopdracht.journal_id_seq'::regclass)"))
    name = Column(String(60), nullable=False)


class Alias(Base):
    __tablename__ = 'alias'
    __table_args__ = {'schema': 'eindopdracht'}

    id = Column(Integer, primary_key=True, server_default=text("nextval('eindopdracht.alias_id_seq'::regclass)"))
    hgnc_symbol = Column(String(15), nullable=False)
    gene_id = Column(ForeignKey('eindopdracht.gene.id'), nullable=False)

    gene = relationship('Gene')


class Article(Base):
    __tablename__ = 'article'
    __table_args__ = {'schema': 'eindopdracht'}

    id = Column(Integer, primary_key=True, server_default=text("nextval('eindopdracht.article_id_seq'::regclass)"))
    title = Column(String(200), nullable=False)
    pubmed_id = Column(Integer, nullable=False)
    doi = Column(String(60), nullable=False)
    publication_date = Column(Date, nullable=False)
    abstract = Column(String(3000), nullable=False)
    journal_id = Column(ForeignKey('eindopdracht.journal.id'), nullable=False)

    journal = relationship('Journal')
    genes = relationship('Gene', secondary='eindopdracht.article_gene')


t_genepanel_gene = Table(
    'genepanel_gene', metadata,
    Column('gene_id', ForeignKey('eindopdracht.gene.id'), nullable=False),
    Column('genepanel_id', ForeignKey('eindopdracht.genepanel.id'), nullable=False),
    schema='eindopdracht'
)


class GenepanelSymbol(Base):
    __tablename__ = 'genepanel_symbol'
    __table_args__ = {'schema': 'eindopdracht'}

    id = Column(Integer, primary_key=True,
                server_default=text("nextval('eindopdracht.genepanel_symbol_id_seq'::regclass)"))
    symbol = Column(String(15), nullable=False)
    gene_id = Column(ForeignKey('eindopdracht.gene.id'), nullable=False)

    gene = relationship('Gene')


t_article_gene = Table(
    'article_gene', metadata,
    Column('gene_id', ForeignKey('eindopdracht.gene.id'), nullable=False),
    Column('article_id', ForeignKey('eindopdracht.article.id'), nullable=False),
    schema='eindopdracht'
)
