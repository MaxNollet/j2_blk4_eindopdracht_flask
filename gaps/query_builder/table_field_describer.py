from gaps.models import *


class TableFieldDescriber:
    """A class which groups methods that specify how fields
       can be selected and from which table the field is
       from.
    """
    column_name_describer = dict()

    def __init__(self):
        self.table_columns_gene()
        self.table_columns_genepanel_symbol()
        self.table_columns_alias()
        self.table_columns_query()
        self.table_columns_article()
        self.table_columns_journal()
        self.table_columns_disease()
        self.table_columns_genepanel()
        self.table_columns_inheritance_type()

    def table_columns_gene(self) -> None:
        """A method which describes how fields from the table
           'gene' can be selected.
        """
        table = Gene.__tablename__
        self.column_name_describer["Gene ID"] = {
            "table": table,
            "column": Gene.id
        }
        self.column_name_describer["NCBI gene ID"] = {
            "table": table,
            "column": Gene.ncbi_gene_id
        }
        self.column_name_describer["Gene symbol"] = {
            "table": table,
            "column": Gene.hgnc_symbol
        }
        self.column_name_describer["In genepanel"] = {
            "table": table,
            "column": Gene.in_genepanel
        }
        return None

    def table_columns_genepanel_symbol(self) -> None:
        """A method which describes how fields from the table
           'genepanel_symbol' can be selected.
        """
        table = GenepanelSymbol.__tablename__
        self.column_name_describer["Genepanel symbol ID"] = {
            "table": table,
            "column": GenepanelSymbol.id
        }
        self.column_name_describer["Genepanel symbol"] = {
            "table": table,
            "column": GenepanelSymbol.symbol
        }
        return None

    def table_columns_alias(self) -> None:
        """A method which describes how fields from the table
           'alias' can be selected.
        """
        table = t_gene_alias.description
        self.column_name_describer["Alias ID"] = {
            "table": table,
            "column": Alias.id
        }
        self.column_name_describer["Alias Symbol"] = {
            "table": table,
            "column": Alias.hgnc_symbol
        }
        return None

    def table_columns_query(self) -> None:
        """A method which describes how fields from the table
           'query' can be selected.
        """
        table = Query.__tablename__
        self.column_name_describer["Query ID"] = {
            "table": table,
            "column": Query.id
        }
        self.column_name_describer["Query"] = {
            "table": table,
            "column": Query.query
        }
        return None

    def table_columns_article(self) -> None:
        """A method which describes how fields from the table
            'article' can be selected.
        """
        table = Article.__tablename__
        self.column_name_describer["Article ID"] = {
            "table": table,
            "column": Article.id
        }
        self.column_name_describer["Article title"] = {
            "table": table,
            "column": Article.title
        }
        self.column_name_describer["Article PubMed ID"] = {
            "table": table,
            "column": Article.pubmed_id
        }
        self.column_name_describer["Article DOI"] = {
            "table": table,
            "column": Article.doi
        }
        self.column_name_describer["Article publication date"] = {
            "table": table,
            "column": Article.publication_date
        }
        self.column_name_describer["Article abstract"] = {
            "table": table,
            "column": Article.abstract
        }
        return None

    def table_columns_journal(self) -> None:
        """A method which describes how fields from the table
           'journal' can be selected.
        """
        table = Journal.__tablename__
        self.column_name_describer["Journal ID"] = {
            "table": table,
            "column": Journal.id
        }
        self.column_name_describer["Journal name"] = {
            "table": table,
            "column": Journal.name
        }
        return None

    def table_columns_disease(self) -> None:
        """A method which describes how fields from the table
           'disease' can be selected.
        """
        table = Disease.__tablename__
        self.column_name_describer["Disease ID"] = {
            "table": table,
            "column": Disease.id
        }
        self.column_name_describer["Disease MESH ID"] = {
            "table": table,
            "column": Disease.mesh_id
        }
        self.column_name_describer["Disease"] = {
            "table": table,
            "column": Disease.disease
        }
        return None

    def table_columns_genepanel(self) -> None:
        """A method which describes how fields from the table
           'genepanel' can be selected.
        """
        table = Genepanel.__tablename__
        self.column_name_describer["Genepanel ID"] = {
            "table": table,
            "column": Genepanel.id
        }
        self.column_name_describer["Genepanel abbreviation"] = {
            "table": table,
            "column": Genepanel.abbreviation
        }
        return None

    def table_columns_inheritance_type(self) -> None:
        """A method which describes how fields from the table
           'inheritance_type' can be selected.
        """
        table = InheritanceType.__tablename__
        self.column_name_describer["Inheritance type ID"] = {
            "table": table,
            "column": InheritanceType.id
        }
        self.column_name_describer["Inheritance type"] = {
            "table": table,
            "column": InheritanceType.type
        }
        return None

    def get_field_names(self):
        """A method which returns a list of all fields that
           can be selected from the database.

        :return All valid selectable fields from different tables (Tuple[str]).
        """
        return sorted(self.column_name_describer.keys())
