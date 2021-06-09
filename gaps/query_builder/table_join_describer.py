from gaps.models import *


class TableJoinDescriber:
    table_joins = dict()

    def __init__(self):
        self.table_genepanel_symbol()
        self.table_alias()
        self.table_query()
        self.table_gene()
        self.table_journal()
        self.table_disease()
        self.table_genepanel()
        self.table_inheritance_type()
        self.table_t_query_article()
        self.table_t_article_gene()
        self.table_t_article_disease()
        self.table_t_gene_alias()
        self.table_t_genepanel_gene()
        self.t_genepanel_inheritance()

    def table_genepanel_symbol(self) -> None:
        self.table_joins[GenepanelSymbol.__tablename__] = {
            "required": Gene.__tablename__,
            "table": GenepanelSymbol,
            "id1": Gene.genepanel_symbol_id,
            "id2": GenepanelSymbol.id
        }
        return None

    def table_alias(self) -> None:
        self.table_joins[Alias.__tablename__] = {
            "required": (t_gene_alias.description,),
            "table": Alias,
            "id1": t_gene_alias.c.alias_id,
            "id2": Alias.id
        }
        return None

    def table_query(self) -> None:
        self.table_joins[Query.__tablename__] = {
            "required": (t_query_article.description,),
            "table": Query,
            "id1": t_query_article.c.query_id,
            "id2": Query.id
        }
        return None

    def table_gene(self) -> None:
        self.table_joins[Gene.__tablename__] = {
            "required": (t_article_gene.description,),
            "table": Gene,
            "id1": t_article_gene.c.gene_id,
            "id2": Gene.id
        }

    def table_journal(self) -> None:
        self.table_joins[Journal.__tablename__] = {
            "table": Journal,
            "id1": Article.journal_id,
            "id2": Journal.id
        }
        return None

    def table_disease(self) -> None:
        self.table_joins[Disease.__tablename__] = {
            "required": (t_article_disease.description,),
            "table": Disease,
            "id1": t_article_disease.c.disease_id,
            "id2": Disease.id
        }
        return None

    def table_genepanel(self) -> None:
        self.table_joins[Genepanel.__tablename__] = {
            "required": (t_genepanel_gene.description,),
            "table": Genepanel,
            "id1": t_genepanel_gene.c.genepanel_id,
            "id2": Genepanel.id
        }
        return None

    def table_inheritance_type(self) -> None:
        self.table_joins[InheritanceType.__tablename__] = {
            "required": (t_genepanel_inheritance.description,),
            "table": InheritanceType,
            "id1": t_genepanel_inheritance.c.inheritance_type_id,
            "id2": InheritanceType.id
        }
        return None

    def table_t_query_article(self) -> None:
        self.table_joins[t_query_article.description] = {
            "table": t_query_article,
            "id1": t_query_article.c.article_id,
            "id2": Article.id
        }
        return None

    def table_t_article_gene(self) -> None:
        self.table_joins[t_article_gene.description] = {
            "table": t_article_gene,
            "id1": t_article_gene.c.article_id,
            "id2": Article.id
        }
        return None

    def table_t_article_disease(self) -> None:
        self.table_joins[t_article_disease.description] = {
            "table": t_article_disease,
            "id1": t_article_disease.c.article_id,
            "id2": Article.id
        }
        return None

    def table_t_gene_alias(self) -> None:
        self.table_joins[t_gene_alias.description] = {
            "required": (Gene.__tablename__,),
            "table": t_gene_alias,
            "id1": t_gene_alias.c.gene_id,
            "id2": Gene.id
        }
        return None

    def table_t_genepanel_gene(self) -> None:
        self.table_joins[t_genepanel_gene.description] = {
            "required": (Gene.__tablename__,),
            "table": t_genepanel_gene,
            "id1": t_genepanel_gene.c.gene_id,
            "id2": Gene.id
        }
        return None

    def t_genepanel_inheritance(self) -> None:
        self.table_joins[t_genepanel_inheritance.description] = {
            "required": (Genepanel.__tablename__,),
            "table": t_genepanel_inheritance.description,
            "id1": t_genepanel_inheritance.c.genepanel_id,
            "id2": Genepanel.id
        }
