# coding: utf-8
from dataclasses import dataclass

from flask_sqlalchemy import SQLAlchemy
from sqlalchemy import Boolean, Column, Date, ForeignKey, Integer, String, Table, Text, text
from sqlalchemy.dialects.postgresql import UUID
from sqlalchemy.orm import relationship

db = SQLAlchemy()

Model = db.Model
metadata = db.metadata


@dataclass
class Alias(Model):
    """A class which maps to the table 'alias'
       in the database.
    """
    __tablename__ = 'alias'
    __table_args__ = {'schema': 'eindopdracht'}

    id = Column(Integer, primary_key=True, server_default=text("nextval('eindopdracht.alias_id_seq'::regclass)"))
    hgnc_symbol = Column(String(30), nullable=False, unique=True)

    genes = relationship('Gene', secondary='eindopdracht.gene_alias')


@dataclass
class Gene(Model):
    """A class which maps to the table 'gene' in
       the database.
    """
    __tablename__ = 'gene'
    __table_args__ = {'schema': 'eindopdracht'}

    id = Column(Integer, primary_key=True, server_default=text("nextval('eindopdracht.gene_id_seq'::regclass)"))
    ncbi_gene_id = Column(Integer)
    hgnc_symbol = Column(String(30), nullable=False, unique=True)
    in_genepanel = Column(Boolean, nullable=False, server_default=text("false"))
    genepanel_symbol_id = Column(ForeignKey('eindopdracht.genepanel_symbol.id'))

    genepanels = relationship('Genepanel', secondary='eindopdracht.genepanel_gene')
    querys = relationship('Query', secondary='eindopdracht.query_gene')


@dataclass
class Genepanel(Model):
    """A class which maps to the table 'genepanel'
       in the database.
    """
    __tablename__ = 'genepanel'
    __table_args__ = {'schema': 'eindopdracht'}

    id = Column(Integer, primary_key=True, server_default=text("nextval('eindopdracht.genepanel_id_seq'::regclass)"))
    abbreviation = Column(String(40), nullable=False, unique=True)

    inheritance_types = relationship('InheritanceType', secondary='eindopdracht.genepanel_inheritance')


@dataclass
class InheritanceType(Model):
    """A class which maps to the table 'inheritance_type'
       in the database.
    """
    __tablename__ = 'inheritance_type'
    __table_args__ = {'schema': 'eindopdracht'}

    id = Column(Integer, primary_key=True,
                server_default=text("nextval('eindopdracht.inheritance_type_id_seq'::regclass)"))
    type = Column(String(15), nullable=False, unique=True)


@dataclass
class Journal(Model):
    """A class which maps to the table 'journal'
       in the database.
    """
    __tablename__ = 'journal'
    __table_args__ = {'schema': 'eindopdracht'}

    id = Column(Integer, primary_key=True, server_default=text("nextval('eindopdracht.journal_id_seq'::regclass)"))
    name = Column(String(60), nullable=False)


@dataclass
class Option(Model):
    """A class which maps to the table 'option'
       in the database.
    """
    __tablename__ = 'options'
    __table_args__ = {'schema': 'eindopdracht'}

    id = Column(Integer, primary_key=True, server_default=text("nextval('eindopdracht.options_id_seq'::regclass)"))
    date_after = Column(Date, nullable=False)
    date_before = Column(Date, nullable=False)


@dataclass
class Symbol(Model):
    """A class which maps to the table 'symbol'
       in the database.
    """
    __tablename__ = 'symbol'
    __table_args__ = {'schema': 'eindopdracht'}

    id = Column(Integer, primary_key=True, server_default=text("nextval('eindopdracht.symbol_id_seq'::regclass)"))
    symbol = Column(String(80), nullable=False, unique=True)


@dataclass
class Article(Model):
    """A class which maps to the table 'article'
       in the database.
    """
    __tablename__ = 'article'
    __table_args__ = {'schema': 'eindopdracht'}

    id = Column(Integer, primary_key=True, server_default=text("nextval('eindopdracht.article_id_seq'::regclass)"))
    title = Column(String(200), nullable=False)
    pubmed_id = Column(Integer)
    doi = Column(String(60), nullable=False, unique=True)
    publication_date = Column(Date)
    abstract = Column(String(3000), nullable=False)
    journal_id = Column(ForeignKey('eindopdracht.journal.id'))

    journal = relationship('Journal')
    genes = relationship('Gene', secondary='eindopdracht.article_gene')


@dataclass
class GenepanelSymbol(Model):
    """A class which maps to the table 'genepanel_symbol'
       in the database.
    """
    __tablename__ = 'genepanel_symbol'
    __table_args__ = {'schema': 'eindopdracht'}

    id = Column(Integer, primary_key=True,
                server_default=text("nextval('eindopdracht.genepanel_symbol_id_seq'::regclass)"))
    symbol = Column(String(30), nullable=False, unique=True)
    # gene_id = Column(ForeignKey('eindopdracht.gene.id'), nullable=False)

    gene = relationship('Gene')


@dataclass
class Query(Model):
    """A class which maps to the table 'query'
       in the database.
    """
    __tablename__ = 'query'
    __table_args__ = {'schema': 'eindopdracht'}

    id = Column(UUID, primary_key=True)
    query = Column(Text, nullable=False)
    options_id = Column(ForeignKey('eindopdracht.options.id'))

    options = relationship('Option')
    symbols = relationship('Symbol', secondary='eindopdracht.query_symbol')


t_gene_alias = Table(
    'gene_alias', metadata,
    Column('gene_id', ForeignKey('eindopdracht.gene.id'), nullable=False),
    Column('alias_id', ForeignKey('eindopdracht.alias.id'), nullable=False),
    schema='eindopdracht'
)

t_genepanel_gene = Table(
    'genepanel_gene', metadata,
    Column('gene_id', ForeignKey('eindopdracht.gene.id'), nullable=False),
    Column('genepanel_id', ForeignKey('eindopdracht.genepanel.id'), nullable=False),
    schema='eindopdracht'
)

t_genepanel_inheritance = Table(
    'genepanel_inheritance', metadata,
    Column('genepanel_id', ForeignKey('eindopdracht.genepanel.id'), nullable=False),
    Column('inheritance_type_id', ForeignKey('eindopdracht.inheritance_type.id'), nullable=False),
    schema='eindopdracht'
)

t_article_gene = Table(
    'article_gene', metadata,
    Column('gene_id', ForeignKey('eindopdracht.gene.id'), nullable=False),
    Column('article_id', ForeignKey('eindopdracht.article.id'), nullable=False),
    schema='eindopdracht'
)

t_query_gene = Table(
    'query_gene', metadata,
    Column('query_id', ForeignKey('eindopdracht.query.id'), nullable=False),
    Column('gene_id', ForeignKey('eindopdracht.gene.id'), nullable=False),
    schema='eindopdracht'
)

t_query_symbol = Table(
    'query_symbol', metadata,
    Column('query_id', ForeignKey('eindopdracht.query.id'), nullable=False),
    Column('symbol_id', ForeignKey('eindopdracht.symbol.id'), nullable=False),
    schema='eindopdracht'
)
