# coding: utf-8
import uuid
from dataclasses import dataclass

from flask_sqlalchemy import SQLAlchemy
from sqlalchemy import Boolean, Column, Date, ForeignKey, Integer, String, Table, Text, text
from sqlalchemy.dialects.postgresql import UUID
from sqlalchemy.orm import relationship

db = SQLAlchemy()

Model = db.Model
metadata = db.metadata


def _unique(session, cls, hashfunc, queryfunc, constructor, arg, kw):
    cache = getattr(session, '_unique_cache', None)
    if cache is None:
        session._unique_cache = cache = {}

    key = (cls, hashfunc(*arg, **kw))
    if key in cache:
        return cache[key]
    else:
        with session.no_autoflush:
            q = session.query(cls)
            q = queryfunc(q, *arg, **kw)
            obj = q.first()
            if not obj:
                obj = constructor(*arg, **kw)
                session.add(obj)
        cache[key] = obj
        return obj


class UniqueMixin(object):
    @classmethod
    def unique_hash(cls, *arg, **kw):
        raise NotImplementedError()

    @classmethod
    def unique_filter(cls, query, *arg, **kw):
        raise NotImplementedError()

    @classmethod
    def as_unique(cls, session, *arg, **kw):
        return _unique(
                    session,
                    cls,
                    cls.unique_hash,
                    cls.unique_filter,
                    cls,
                    arg, kw
               )


@dataclass
class Gene(UniqueMixin, Model):
    """A class which maps to the table 'gene' in
       the database.
    """
    __tablename__ = 'gene'
    __table_args__ = {'schema': 'eindopdracht'}

    id: int = Column(Integer, primary_key=True, server_default=text("nextval('eindopdracht.gene_id_seq'::regclass)"))
    ncbi_gene_id: int = Column(Integer)
    hgnc_symbol: str = Column(String(30), nullable=False, unique=True)
    in_genepanel: bool = Column(Boolean, nullable=False, server_default=text("false"))

    queries: list = relationship('Query', secondary='eindopdracht.query_gene')
    genepanels: list = relationship('Genepanel', secondary='eindopdracht.genepanel_gene')
    aliases: list = relationship('Alias', secondary='eindopdracht.gene_alias', back_populates='genes')

    @classmethod
    def unique_hash(cls, hgnc_symbol):
        return hgnc_symbol

    @classmethod
    def unique_filter(cls, query, hgnc_symbol):
        return query.filter(Gene.hgnc_symbol == hgnc_symbol)


@dataclass
class Genepanel(Model):
    """A class which maps to the table 'genepanel'
       in the database.
    """
    __tablename__ = 'genepanel'
    __table_args__ = {'schema': 'eindopdracht'}

    id: int = Column(Integer, primary_key=True,
                     server_default=text("nextval('eindopdracht.genepanel_id_seq'::regclass)"))
    abbreviation: str = Column(String(40), nullable=False, unique=True)

    inheritance_types: list = relationship('InheritanceType', secondary='eindopdracht.genepanel_inheritance')


@dataclass
class InheritanceType(Model):
    """A class which maps to the table 'inheritance_type'
       in the database.
    """
    __tablename__ = 'inheritance_type'
    __table_args__ = {'schema': 'eindopdracht'}

    id: int = Column(Integer, primary_key=True,
                     server_default=text("nextval('eindopdracht.inheritance_type_id_seq'::regclass)"))
    type: str = Column(String(15), nullable=False, unique=True)


@dataclass
class Journal(Model):
    """A class which maps to the table 'journal'
       in the database.
    """
    __tablename__ = 'journal'
    __table_args__ = {'schema': 'eindopdracht'}

    id: int = Column(Integer, primary_key=True, server_default=text("nextval('eindopdracht.journal_id_seq'::regclass)"))
    name: str = Column(String(60), nullable=False)


@dataclass
class Option(Model):
    """A class which maps to the table 'option'
       in the database.
    """
    __tablename__ = 'options'
    __table_args__ = {'schema': 'eindopdracht'}

    id: int = Column(Integer, primary_key=True, server_default=text("nextval('eindopdracht.options_id_seq'::regclass)"))
    date_after: Date = Column(Date, nullable=False)
    date_before: Date = Column(Date, nullable=False)


@dataclass
class Symbol(Model):
    """A class which maps to the table 'symbol'
       in the database.
    """
    __tablename__ = 'symbol'
    __table_args__ = {'schema': 'eindopdracht'}

    id: int = Column(Integer, primary_key=True, server_default=text("nextval('eindopdracht.symbol_id_seq'::regclass)"))
    symbol: str = Column(String(80), nullable=False, unique=True)


@dataclass
class Alias(UniqueMixin, Model):
    """A class which maps to the table 'alias'
       in the database.
    """
    __tablename__ = 'alias'
    __table_args__ = {'schema': 'eindopdracht'}

    id: int = Column(Integer, primary_key=True, server_default=text("nextval('eindopdracht.alias_id_seq'::regclass)"))
    hgnc_symbol: str = Column(String(30), nullable=False, unique=True)

    genes: list = relationship('Gene', secondary='eindopdracht.gene_alias', back_populates='aliases')

    @classmethod
    def unique_hash(cls, hgnc_symbol):
        return hgnc_symbol

    @classmethod
    def unique_filter(cls, query, hgnc_symbol):
        return query.filter(Alias.hgnc_symbol == hgnc_symbol)


@dataclass
class Article(Model):
    """A class which maps to the table 'article'
       in the database.
    """
    __tablename__ = 'article'
    __table_args__ = {'schema': 'eindopdracht'}

    id: int = Column(Integer, primary_key=True, server_default=text("nextval('eindopdracht.article_id_seq'::regclass)"))
    title: str = Column(String(200), nullable=False)
    pubmed_id: int = Column(Integer)
    doi: str = Column(String(60), nullable=False, unique=True)
    publication_date: Date = Column(Date)
    abstract: str = Column(String(3000), nullable=False)
    journal_id: int = Column(ForeignKey('eindopdracht.journal.id'))

    journal: Journal = relationship('Journal')
    genes: list = relationship('Gene', secondary='eindopdracht.article_gene')


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


@dataclass
class GenepanelSymbol(Model):
    """A class which maps to the table 'genepanel_symbol'
       in the database.
    """
    __tablename__ = 'genepanel_symbol'
    __table_args__ = {'schema': 'eindopdracht'}

    id: int = Column(Integer, primary_key=True,
                     server_default=text("nextval('eindopdracht.genepanel_symbol_id_seq'::regclass)"))
    symbol: str = Column(String(30), nullable=False, unique=True)
    gene_id: int = Column(ForeignKey('eindopdracht.gene.id'), nullable=False)

    gene: Gene = relationship('Gene')


@dataclass
class Query(Model):
    """A class which maps to the table 'query'
       in the database.
    """
    __tablename__ = 'query'
    __table_args__ = {'schema': 'eindopdracht'}

    id: uuid.uuid4 = Column(UUID, primary_key=True, default=uuid.uuid4)
    query: str = Column(Text, nullable=False)
    options_id: int = Column(ForeignKey('eindopdracht.options.id'))

    options: Option = relationship('Option')
    symbols: list = relationship('Symbol', secondary='eindopdracht.query_symbol')


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
