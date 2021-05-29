from sqlalchemy import select, bindparam

from gaps.models import Gene, Alias, GenepanelSymbol, Genepanel, InheritanceType


class SelectStatements:
    """A class which groups together several statements
       for selecting values from different tables in the
       database.
    """

    @staticmethod
    def _select_gene():
        """A statement which selects hgnc_symbols and ids
           from the gene-table only if hgnc_symbol is
           in the filter.

        :return Select-statement for specific values from the gene-table (Select).
        """
        return select(
            Gene.hgnc_symbol, Gene.id
        ).where(
            Gene.hgnc_symbol.in_(bindparam("values"))
        )

    @staticmethod
    def _select_alias():
        """A statement wich selects hgnc_symbols and ids
           from the alias-table only if hgnc_symbol is
           in the filter.

        :return Select-statement for specific values from the alias-table (Select).
        """
        return select(
            Alias.hgnc_symbol, Alias.id
        ).where(
            Alias.hgnc_symbol.in_(bindparam("values"))
        )

    @staticmethod
    def _select_genepanel_symbol():
        """A statement which selects symbols and ids
           from the genepanel-table only if symbol is
           in the filter.

        :return Select-statement for specific values from the genepanel-table (Select).
        """
        return select(
            GenepanelSymbol.symbol, GenepanelSymbol.id
        ).where(
            GenepanelSymbol.symbol.in_(bindparam("values"))
        )

    @staticmethod
    def _select_genepanel():
        """A statement which selects abbreviations and ids
           from the genepanel-table only if abbreviation is
           in the filter.

        :return Select-statement for specific values from the genepanel-table (Select).
        """
        return select(
            Genepanel.abbreviation, Genepanel.id
        ).where(
            Genepanel.abbreviation.in_(bindparam("values"))
        )

    @staticmethod
    def _select_inheritance_type():
        """A statement which selects types and ids from
           the inheritance-table only when type is
           in the filter.

        :return Select-statement for specific values form the inheritance-table (Select).
        """
        return select(
            InheritanceType.type, InheritanceType.id
        ).where(
            InheritanceType.type.in_(bindparam("values"))
        )
