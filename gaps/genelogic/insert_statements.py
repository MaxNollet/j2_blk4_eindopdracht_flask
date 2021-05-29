from sqlalchemy.dialects.postgresql import insert

from gaps.models import Gene, Alias, GenepanelSymbol, Genepanel, InheritanceType, t_gene_alias, t_genepanel_gene, \
    t_genepanel_inheritance


class InsertStatements:
    """A class which groups together several statements
       for inserting values into different tables into
       the database.
    """

    @staticmethod
    def _insert_gene():
        """A statement for inserting values into the
           gene-table.

        :return Insert-statement for the gene-table (Insert).
        """
        return insert(Gene).on_conflict_do_nothing()

    @staticmethod
    def _insert_alias():
        """A statement for inserting values into the
           alias-table.

        :return Insert-statement for the alias-table (Insert).
        """
        return insert(Alias).on_conflict_do_nothing()

    @staticmethod
    def _insert_genepanel_symbol():
        """A statement for inserting values into the
           genepanel_symbol-table.

        :return Insert statement for the genepanel_symbol-table (Insert).
        """
        return insert(GenepanelSymbol).on_conflict_do_nothing()

    @staticmethod
    def _insert_genepanel():
        """A statement for inserting values into the
           genepanel-table.

        :return Insert-statement for the genepanel-table (Insert).
        """
        return insert(Genepanel).on_conflict_do_nothing()

    @staticmethod
    def _insert_inheritance_type():
        """A statement for inserting values into the
           inheritance_type-table.

        :return Insert-statement for the inheritance_type-table (Insert).
        """
        return insert(InheritanceType).on_conflict_do_nothing()

    @staticmethod
    def _insert_relation_gene_alias():
        """A statement for inserting values into the
           gene_alias-table.

        :return Insert-statement for the gene_alias-table (Insert).
        """
        return insert(t_gene_alias).on_conflict_do_nothing()

    @staticmethod
    def _insert_relation_gene_genepanel():
        """A statement for inserting values into the
           genepanel_gene-table.

        :return Insert-statement for the genepanel_gene-table (Insert).
        """
        return insert(t_genepanel_gene).on_conflict_do_nothing()

    @staticmethod
    def _insert_relation_genepanel_inheritance():
        """A statement for inserting values into the
           genepanel_inheritance-table.

        :return Insert-statement for the genepanel_inheritance-table (Insert).
        """
        return insert(t_genepanel_inheritance).on_conflict_do_nothing()
