from sqlalchemy.dialects.postgresql import insert

from gaps.genelogic.statement_groups import statement_group
from gaps.models import Gene, Alias, GenepanelSymbol, Genepanel, \
    InheritanceType, t_gene_alias, t_genepanel_gene, \
    t_genepanel_inheritance, Article, Journal


class InsertStatements:
    """A class which groups together several statements
       for inserting values into different tables into
       the database.
    """

    @staticmethod
    @statement_group(table="gene")
    def _insert_gene():
        """A statement for inserting values into the
           gene-table.

        :return Insert-statement for the gene-table (Insert).
        """
        return insert(Gene).on_conflict_do_nothing()

    @staticmethod
    @statement_group(table="alias")
    def _insert_alias():
        """A statement for inserting values into the
           alias-table.

        :return Insert-statement for the alias-table (Insert).
        """
        return insert(Alias).on_conflict_do_nothing()

    @staticmethod
    @statement_group(table="genepanel_symbol")
    def _insert_genepanel_symbol():
        """A statement for inserting values into the
           genepanel_symbol-table.

        :return Insert statement for the genepanel_symbol-table (Insert).
        """
        return insert(GenepanelSymbol).on_conflict_do_nothing()

    @staticmethod
    @statement_group(table="genepanel")
    def _insert_genepanel():
        """A statement for inserting values into the
           genepanel-table.

        :return Insert-statement for the genepanel-table (Insert).
        """
        return insert(Genepanel).on_conflict_do_nothing()

    @staticmethod
    @statement_group(table="inheritance_type")
    def _insert_inheritance_type():
        """A statement for inserting values into the
           inheritance_type-table.

        :return Insert-statement for the inheritance_type-table (Insert).
        """
        return insert(InheritanceType).on_conflict_do_nothing()

    @staticmethod
    @statement_group(table="gene_alias")
    def _insert_relation_gene_alias():
        """A statement for inserting values into the
           gene_alias-table.

        :return Insert-statement for the gene_alias-table (Insert).
        """
        return insert(t_gene_alias).on_conflict_do_nothing()

    @staticmethod
    @statement_group(table="genepanel_gene")
    def _insert_relation_gene_genepanel():
        """A statement for inserting values into the
           genepanel_gene-table.

        :return Insert-statement for the genepanel_gene-table (Insert).
        """
        return insert(t_genepanel_gene).on_conflict_do_nothing()

    @staticmethod
    @statement_group(table="genepanel_inheritance")
    def _insert_relation_genepanel_inheritance():
        """A statement for inserting values into the
           genepanel_inheritance-table.

        :return Insert-statement for the genepanel_inheritance-table (Insert).
        """
        return insert(t_genepanel_inheritance).on_conflict_do_nothing()

    @staticmethod
    @statement_group(table="article")
    def _insert_article():
        return insert(Article).on_conflict_do_nothing()
