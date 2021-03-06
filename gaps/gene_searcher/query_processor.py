# import ssl
import uuid
from dataclasses import dataclass, field
from datetime import datetime
from os import environ
from typing import Dict
from xml.etree import ElementTree as etree
import requests
from Bio import Entrez
from gaps.genelogic.database_inserter import DatabaseInserter

Entrez.email = environ.get("EMAIL_ENTREZ")
Entrez.email = "mjh.nollet@student.han.nl"


# def main():
#     # ssl._create_default_https_context = ssl._create_unverified_context
#
#     query = "((variant [tiab] OR variants [tiab] OR mutation [tiab] OR mutations [tiab] OR substitutions [tiab] OR substitution [tiab] ) \
#     AND (\"loss of function\" [tiab] OR \"loss-of-function\" [tiab] OR \"haplo-insufficiency\" [tiab] OR haploinsufficiency [tiab] \
#     OR \"bi-allelic\" [tiab] OR \"biallelic\" [tiab] OR recessive [tiab] OR homozygous [tiab] OR heterozygous [tiab] OR \"de novo\" \
#     [tiab] OR dominant [tiab] OR \" X-linked\" [tiab]) AND (\"intellectual\" [tiab] OR \"mental retardation\" [tiab] OR \"cognitive\" \
#     [tiab] OR \"developmental\" [tiab] OR \"neurodevelopmental\" [tiab]) AND “last 2 years”[dp] AND KDM3B)"
#
#     print(Entrez.email)
#     query = "((\"2021\"[Date - Publication] : \"3000\"[Date - Publication])) AND (CDH8[Text Word])"
#     # ((Loss of sight[ALL]) AND (Blindness[ALL])) AND (Loss of hearing[ALL])


class GeneSearcher:
    def __init__(self):
        self.search_results = None

        # disabled
        # =-=-=-=-=-=
        self.specify_gene = None  # true or false to use it or not
        self.include_exclude = None  # true or false
        self.specific_gene_symbols = None  # list to include or exclude
        # =-=-=-=-=-=
        # creates a new uuid for each new search
        self.uuid_query = uuid.uuid4()

        self._cached_values: Dict[str, set] = dict()
        # checks if a gene is already in the cache, but still
        # inserts in the between tables
        tables = ["gene", "disease"]
        for table in tables:
            self._cached_values[table] = set()
        self.db = InsertDB()

    def fetch_results(self, parameters: dict) -> int:
        """A method which queries PubMed and returns the count
           of results matching the search parameters. Can be used
           to quickly inform the user if their inserted parameters
           have any results of not.

        :param parameters Parameters entered by the user (dict).
        :return Count of matching articles (int).
        """  # query from querybuilder
        query = parameters.get("input_generated_query")
        min_date: datetime = parameters.get("input_date_after")
        max_date: datetime = parameters.get("input_date_before")
        self.query_validator(query)
        self.db.query_list.append({"id": str(self.uuid_query), "query": query})
        # options_id moet nog verwerkt worden TODO
        if any(date is not None for date in (min_date, max_date)):
            if min_date is None:  # if there isn't a min date
                raise NoDateAfterSpecified()
            if max_date is None:
                raise NoDateBeforeSpecified()
            self.db.query_options_list.append(  # structure for db
                {"date_before": min_date, "date_after": max_date})
            search_results = Entrez.read(
                Entrez.esearch(
                    db="pubmed",
                    term=query,  # search pubmed with query
                    usehistory="y",
                    mindate=min_date.strftime("%Y/%m/%d"),
                    maxdate=max_date.strftime("%Y/%m/%d")
                )
            )
        else:
            search_results = Entrez.read(
                Entrez.esearch(
                    db="pubmed",
                    term=query,
                    usehistory="y"  # the Entrez history feature
                )
            )
        self.search_results = search_results
        return int(search_results["Count"])

    def results_query(self):
        """
        Search for genes and diseases in pubmed and insert the found data
        in the database.
        :return:
        """
        self.query_pubmed()  # looks up data from pubmed with the query
        # print(self.db, "self.db")  # example ouput for in db
        # print(self.db.disease_list, "disease list")
        # print(self.db.genes_list, "genes_list")
        # print(self.db.article_gene, "article gene")
        if not self.db.genes_list:  # no genes found.
            raise NoGeneFound
        else:  # found gene and inserts into database
            db = DatabaseInserter()
            db.insert_search_results(self.db)

    @staticmethod
    def query_validator(query):
        """
        Validates the query on parentheses being closed off correctly.
        Also checks if search terms contain quotation marks.
        :param query: String which contains the query used as search term
        :return:
        """
        open_list = ["[", "("]
        close_list = ["]", ")"]
        stack = []
        if query is not None or query == "":
            for i in query:
                if i in open_list:
                    stack.append(i)
                elif i in close_list:
                    pos = close_list.index(i)
                    if ((len(stack) > 0) and
                            (open_list[pos] == stack[len(stack) - 1])):
                        stack.pop()
                    else:
                        raise MalformedQuery
            if len(stack) != 0:
                raise MalformedQuery
        else:
            raise NoQuerySpecified

    def query_pubmed(self):
        """
        Uses the query to look for article on pubmed and sets the
        needed information like Title, PMID, doi, etc. And sends the
        pubmed id to pubtator for gene information
        """

        count = int(self.search_results["Count"])
        print("Found %i results" % count)
        batch_size = 1000
        for start in range(0, count, batch_size):
            handle = Entrez.efetch(
                db="pubmed",
                rettype="medline",
                retmode="xml",
                retstart=start,
                retmax=batch_size,
                webenv=self.search_results["WebEnv"],
                query_key=self.search_results["QueryKey"],
            )
            records = Entrez.read(handle)
            for record in records["PubmedArticle"]:
                pubmed_id = record.get("MedlineCitation").get("PMID")
                title = record.get("MedlineCitation").get("Article").get(
                    "ArticleTitle")
                doi_element = record.get("MedlineCitation").get("Article").get(
                    "ELocationID")
                # [StringElement('21479269', attributes={'IdType':
                # 'pubmed'}), StringElement('10.1371/journal.
                # pone.0015669', attributes={'IdType': 'doi'}),
                # StringElement('PMC3066203', attributes={'IdType':
                # 'pmc'})]
                if not doi_element:  # if there isn't a doi number
                    doi = str(uuid.uuid4())
                    # doi = record.get("PubmedData").get("ArticleIdList")[-1]
                    # if "/" not in doi and "." not in doi:
                    #     raise IncorrectArticleFound
                else:
                    if doi_element is not None:
                        try:
                            doi = doi_element[0]
                        except IndexError:
                            doi = None
                publication_date = self.extract_date(record)
                abstract_element = \
                    record.get("MedlineCitation").get("Article").get(
                        "Abstract")
                abstract = None
                if abstract_element is not None:
                    abstract = abstract_element.get("AbstractText")[0]
                journal_name = \
                    record.get("MedlineCitation").get("Article").get(
                        "Journal").get("Title")
                self.db.article_list.append({"title": title,
                                             "pubmed_id": pubmed_id,
                                             "doi": doi, "publication_date":
                                                 publication_date, "abstract":
                                                 abstract, "journal_id":
                                                 journal_name})
                self.db.journal_list.append({"name": journal_name})
                self.db.journal_pk_list.append({  # structure to insert
                    "id": journal_name})  # into the database
        self.pubtator_output()

    @staticmethod
    def extract_date(record: dict):
        """A method which safely extracts a publication date
           from an article, in the correct format for the database.
        :param record Article possibly containing a publication date (Dict).
        :return Extracted date from the article if available.
        """
        pub_date = record.get("MedlineCitation").get("Article") \
            .get("Journal").get("JournalIssue").get("PubDate")
        # DictElement({'Volume': '105', 'Issue': '1-2', 'PubDate': {
        # 'MedlineDate': '1999 Jul-Aug'}}, attributes={'CitedMedium': 'Print'})
        if pub_date is not None:
            year: str = pub_date.get("Year")
            month: str = pub_date.get("Month")
            day: str = pub_date.get("Day")
            date = None
            if year and month and day:
                concatenated = f"{year}-{month}-{day}"
                if month.isalpha():
                    date = datetime.strptime(concatenated, "%Y-%b-%d")
                else:
                    date = datetime.strptime(concatenated, "%Y-%m-%d")
            elif year and month:
                concatenated = f"{year}-{month}"
                if month.isalpha():
                    date = datetime.strptime(concatenated, "%Y-%b")
                else:
                    date = datetime.strptime(concatenated, "%Y-%m")
            elif year:
                date = datetime.strptime(year, "%Y")
            else:
                try:
                    medline_date = pub_date.get("MedlineDate")
                    if medline_date:
                        date = datetime.strptime(medline_date.split("-")[0],
                                                 "%Y %b")
                except ValueError:
                    date = None
            return date
        return None

    @staticmethod
    def url_maker(idlist: tuple):
        """
        Creates the URL used for looking up genes on Pubtator
        :param idlist: List that contains all of the unique ids found in
        query_pubmed
        :return: complete_url: url for Pubtator
        """
        # post request 1000
        # get request = 100
        url = ""
        for pmid in idlist:  # id = pubmedid for article pubtator
            if pmid != idlist[len(idlist) - 1]:
                url += pmid + ","
            else:
                url += pmid
        # url = ",".join(idlist)
        complete_url = "https://www.ncbi.nlm.nih.gov/research/pubtator-api" \
                       "/publications/export/biocxml?pmids={}&" \
                       "concepts=gene,disease".format(url)

        # complete_url is url for pubtator from pubtator API
        return complete_url

    def pubtator_output(self):
        """
        Uses the URL made in url_maker to look up genes on pubtator.
        It contains the genes found in the articles
        :return: articles updated with genes from pubtator
        """
        batch_size = 1000
        count_articles = len(self.db.article_list)
        if count_articles > batch_size:
            for slice_start in range(0, count_articles, batch_size):
                if slice_start >= count_articles - batch_size:
                    slice_end = count_articles
                else:
                    slice_end = slice_start + batch_size - 1
                batch = self.db.article_list[slice_start:slice_end]
                self.pubtator_results_processor(batch)
        else:
            self.pubtator_results_processor(self.db.article_list)

    def pubtator_results_processor(self, batch):
        """
        Get results from pubtator by a post request, and calls function
        to insert in lists from the dataclass 'InsertDB'
        :param batch: size of the request (max 1000 per request in post)
        :return:
        """
        idlist = {}  # dict with pubmed_id and doi (or uuid if no doi)
        for article in batch:
            idlist[article["pubmed_id"]] = str(article["doi"])
        json_parameters = {"pmids": [str(element)
                                     for element in idlist.keys()]}
        json_parameters["concepts"] = ("gene", "disease")
        base_url = "https://www.ncbi.nlm.nih.gov/research/pubtator-api" \
                   "/publications/export/"
        result = requests.post(base_url + "biocxml",
                               json=json_parameters)  # get xml-page pubtator
        if result.status_code == 200:
            print("Request succesful.")
            data = self.parse_results(result)
            for article in idlist.keys():  # article is a int pb_id
                # 0 is always gene, 1 is always diseases
                try:
                    for id_gene, gene in data[article][0].items():
                        if ";" not in id_gene:  # doesn't work with db
                            # if self.specify_gene:  # true to activate
                            #     self._check_include_exclude(id_gene, gene,
                            #                                 article)
                            # else:  # don't use it
                            self._extract_gene(id_gene, gene, article)
                            self.db.query_article.append(  # TEST
                                {"query_id": self.db.query_list[0][
                                    "query"],
                                 "article_id": int(article)})
                    for id_disease, disease in data[article][1].items():
                        self._extract_disease(id_disease, disease, article)
                except KeyError as e:
                    print("An article had no genes or diseases ->", e)
        else:
            print("Request not succesful.")

    # def _check_include_exclude(self, id_gene, gene, article):
    #   # developed to include and exlude the specify genes
    #   # but gives problems for the db with the query_article table
    #     if self.include_exclude:  # true include de genen
    #         print(self.specific_gene_symbols)
    #         if gene in self.specific_gene_symbols:
    #             print("include", gene)
    #             self._extract_gene(id_gene, gene, article)
    #
    #     else:
    #         if gene not in self.specific_gene_symbols:
    #             print("exclude", gene)
    #             self._extract_gene(id_gene, gene, article)

    def _extract_gene(self, id_gene, gene, article):
        """A method which checks if a gene is already in the database,
        to prevent on conflict_do_update problemds. But fills the between
        table article_gene
        :param id_gene: id from gene
        :param gene: gene
        :param article: Id from article so pubmed_id
        :return:
        """
        cache = self._cached_values["gene"]
        if gene not in cache:  # check if gene is already in cache
            cache.add(gene)  # adds to set
            self.db.genes_list.append({"ncbi_gene_id":
                                           int(id_gene),
                                       "hgnc_symbol": str(
                                           gene),
                                       "in_genepanel": False})
        # self.db.article_gene.append(
        #     {"article_id": idlist[article],
        #      "gene_id": gene})
        self.db.article_gene.append(
            {"article_id": int(article),
             "gene_id": gene})

    def _extract_disease(self, id_disease, disease, article):
        """A method which checks if a disease is already in the datsbase,
        to prevent on_conflict_do_update problems. But fills the between
        table article_disease
        :param id_disease: MESH_ID
        :param disease: Disease
        :param article: Id from article so pubmed_id
        :return:
        """
        cache = self._cached_values["disease"]
        if disease not in cache:
            cache.add(disease)
            if id_disease[:4] == "OMIM":
                self.db.disease_list.append(
                    # OMIM in disease_id
                    {"mesh_id": id_disease, "disease": disease})
            else:
                if "_" in id_disease:  # to prevent unwanted slicing
                    self.db.disease_list.append(
                        # eg id_disease MESH:D064752_10001
                        {"mesh_id": id_disease[5:-5],  # eg D064752
                         "disease": disease})
                else:
                    self.db.disease_list.append(
                        {"mesh_id": id_disease[5:],
                         "disease": disease})
        # self.db.article_disease.append(
        # {"disease_id": disease,
        # "article_id": idlist[article]})
        self.db.article_disease.append(
            {"disease_id": disease,
             "article_id": int(article)})

    @staticmethod
    def anno_document(documents):
        """
        Adds the ncbi gene id and the gene per article to list, so it can
        later be added to the Article object.
        :param documents: <passage>info</passage> all the info from pubtator
        :return:
        """
        data_document = []
        id_gene = 0  # if there is no ID's
        for document in documents.findall("passage"):
            for doc in document.findall("annotation"):
                for anno in doc:
                    if anno.text is None:
                        data_document.append([anno.attrib])
                    elif anno.text == "None":  # when there is no MESH:
                        data_document.append(
                            [anno.attrib, "MESH:" + str(id_gene)])
                        id_gene += 1  # if there is no ID
                    else:
                        data_document.append([anno.attrib, anno.text])
        return data_document

    @staticmethod
    def article_id(documents):
        """
        Find al articles in the output from pubtator
        :param documents:
        :return:
        """
        for document in documents.findall("id"):
            return document.text

    def parse_results(self, result):
        """
        Parses all the necessary results out of the Pubtator output.
        :param result: xml form from Pubtator output
        :return: data_pubtator:
        list with dict data_pubtator[pmid] = [genes, mesh]
        """
        tree = etree.fromstring(result.text)
        data_doc = {}
        for documents in tree.findall("document"):
            data_doc[self.article_id(documents)] = self.anno_document(
                documents)
        data_pubtator = {}
        for pmid, data in data_doc.items():  # pubmed id, data pubtator
            genes = {}
            mesh = {}
            count = 1000
            for gene in data:  # gene = {'key': 'identifier'}, '5362']
                if count <= 9999:  # unique id voor mesh/ omim
                    try:               # prevents aliases/ synonyms
                        if gene[0]["key"] == "identifier":
                            if gene[1][:4] == "MESH":
                                gene_id = gene[1] + "_" + str(count)
                                mesh[gene_id] = ""
                                count += 1
                            elif gene[1][:4] == "OMIM":
                                gene_id = gene[1] + "_" + str(count)
                                # print(gene_id)  # verification
                                mesh[gene_id] == ""
                                count += 1
                            else:
                                gene_id = gene[1]
                                genes[gene_id] = ""
                        if gene[0]["key"] == "Identifier":  # mesh
                            # if gene[1][:4] == "MESH":
                            gene_id = gene[1]
                            mesh[gene_id] = ""
                    except KeyError:
                        pass
                else:
                    print("wow, that are a lot of unexpected diseases")
                    pass
                if len(gene) >= 2:
                    if not gene[0]:
                        if gene_id[:4] == "MESH" or gene_id[:4] == "OMIM":
                            if gene[1].strip() != "":
                                mesh[gene_id.strip()] = gene[1].strip()
                        else:
                            genes[gene_id] = gene[1]
            data_pubtator[pmid] = [genes, mesh]
            # print(data_pubtator)
        return data_pubtator


@dataclass
class InsertDB:  # class to insert all results into the database
    genes_list: list = field(default_factory=list)

    article_list: list = field(default_factory=list)
    article_gene: list = field(default_factory=list)
    journal_list: list = field(default_factory=list)
    journal_pk_list: list = field(default_factory=list)

    disease_list: list = field(default_factory=list)
    article_disease: list = field(default_factory=list)

    query_list: list = field(default_factory=list)
    query_options_list: list = field(default_factory=list)
    query_article: list = field(default_factory=list)


class MalformedQuery(Exception):
    """Excepton for when a custom query is not correct"""

    def __init__(self):
        super().__init__("Not a correct query! Check your query for "
                         "opening and closing parenthesis, or try to "
                         "create a new query with the querybuilder.")


class NoGeneFound(Exception):
    """If there isn't a gene found with the query"""

    def __init__(self):
        super().__init__("No genes found with this query.")


class IncorrectArticleFound(Exception):
    """If an article don't have a DOI number, to prevent crash database"""

    def __init__(self):
        super().__init__("Found an article with no DOI.")


class NoQuerySpecified(Exception):
    """Exception for when no query is specified."""

    def __init__(self):
        super().__init__("No query specified! Can't perform a search when "
                         "no query is specified.")


class NoDateAfterSpecified(Exception):
    """Exception for when parameter 'date after' is not specified
       but 'date before' is specified.
    """

    def __init__(self):
        super().__init__("No 'date after' specified! 'Date after'-parameter "
                         "is required when 'date before'-parameter is "
                         "specified.")


class NoDateBeforeSpecified(Exception):
    """Exception for when parameter 'date before' is not specified
       but 'date after' is specified.
    """

    def __init__(self):
        super().__init__("No 'date before' specified! 'Date before'-parameter "
                         "is required when 'date after'-parameter is "
                         "specified.")


if __name__ == "__main__":
    main()
