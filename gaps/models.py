# coding: utf-8
from dataclasses import dataclass

from flask_sqlalchemy import SQLAlchemy
from sqlalchemy import Boolean, Column, Date, ForeignKey, Integer, String, \
    Table, text
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship

db = SQLAlchemy()

Base = declarative_base()
metadata = Base.metadata


@dataclass
class Gene(Base):
    """A class which maps to the table 'gene' in
       the database.
    """
    __tablename__ = 'gene'
    __table_args__ = {'schema': 'eindopdracht'}

    id: int = Column(Integer, primary_key=True, server_default=text(
        "nextval('eindopdracht.gene_id_seq'::regclass)"))
    ncbi_gene_id: str = Column(String(60), nullable=False)
    hgnc_symbol: str = Column(String(15), nullable=False)
    in_genepanel: bool = Column(Boolean, nullable=False)

    genepanels = relationship('Genepanel',
                              secondary='eindopdracht.genepanel_gene')


@dataclass
class Genepanel(Base):
    """A class which maps to the table 'genepanel'
       in the database.
    """
    __tablename__ = 'genepanel'
    __table_args__ = {'schema': 'eindopdracht'}

    id: int = Column(Integer, primary_key=True,
                     server_default=text(
                         "nextval('eindopdracht.genepanel_id_seq'::regclass)"))
    afkorting: str = Column(String(40), nullable=False)
    naam: str = Column(String(100), nullable=False)


@dataclass
class Journal(Base):
    """A class which maps to the table 'journal'
       in the database.
    """
    __tablename__ = 'journal'
    __table_args__ = {'schema': 'eindopdracht'}

    id: int = Column(Integer, primary_key=True, server_default=text(
        "nextval('eindopdracht.journal_id_seq'::regclass)"))
    name: str = Column(String(60), nullable=False)


@dataclass
class Alias(Base):
    """A class which maps to the table 'alias'
       in the database.
    """
    __tablename__ = 'alias'
    __table_args__ = {'schema': 'eindopdracht'}

    id: int = Column(Integer, primary_key=True, server_default=text(
        "nextval('eindopdracht.alias_id_seq'::regclass)"))
    hgnc_symbol: str = Column(String(15), nullable=False)
    gene_id: int = Column(ForeignKey('eindopdracht.gene.id'), nullable=False)

    gene = relationship('Gene')


@dataclass
class Article(Base):
    """A class which maps to the table 'article'
       in the database.
    """
    __tablename__ = 'article'
    __table_args__ = {'schema': 'eindopdracht'}

    id: int = Column(Integer, primary_key=True, server_default=text(
        "nextval('eindopdracht.article_id_seq'::regclass)"))
    title: str = Column(String(200), nullable=False)
    pubmed_id: int = Column(Integer, nullable=False)
    doi: str = Column(String(60), nullable=False)
    publication_date: Date = Column(Date, nullable=False)
    abstract: str = Column(String(3000), nullable=False)
    journal_id: int = Column(ForeignKey('eindopdracht.journal.id'),
                             nullable=False)

    journal = relationship('Journal')
    genes = relationship('Gene', secondary='eindopdracht.article_gene')


t_genepanel_gene = Table(
    'genepanel_gene', metadata,
    Column('gene_id', ForeignKey('eindopdracht.gene.id'), nullable=False),
    Column('genepanel_id', ForeignKey('eindopdracht.genepanel.id'),
           nullable=False),
    schema='eindopdracht'
)


@dataclass
class GenepanelSymbol(Base):
    """A class which maps to the table 'genepanel_symbol'
       in the database.
    """
    __tablename__ = 'genepanel_symbol'
    __table_args__ = {'schema': 'eindopdracht'}

    id: int = Column(Integer, primary_key=True,
                     server_default=text(
                         "nextval('eindopdracht.genepanel_symbol_id_seq'::regclass)"))
    symbol: str = Column(String(15), nullable=False)
    gene_id: int = Column(ForeignKey('eindopdracht.gene.id'), nullable=False)

    gene = relationship('Gene')


t_article_gene = Table(
    'article_gene', metadata,
    Column('gene_id', ForeignKey('eindopdracht.gene.id'), nullable=False),
    Column('article_id', ForeignKey('eindopdracht.article.id'),
           nullable=False),
    schema='eindopdracht'
)


@dataclass
class FileInfo:
    gene: Gene
    # panel: Genepanel
    alias: list
    # p_symbol: GenepanelSymbol
