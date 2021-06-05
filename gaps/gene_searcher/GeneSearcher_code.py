import ssl
from dataclasses import dataclass, field
from datetime import datetime
from os import environ
from xml.etree import ElementTree as etree

import requests
from Bio import Entrez

from gaps.genelogic.database_inserter import DatabaseInserter
from gaps.models import Article, Journal

Entrez.email = environ.get("EMAIL_ENTREZ")
Entrez.email = "mjh.nollet@student.han.nl"


def main():
    ssl._create_default_https_context = ssl._create_unverified_context

    query = "((variant [tiab] OR variants [tiab] OR mutation [tiab] OR mutations [tiab] OR substitutions [tiab] OR substitution [tiab] ) \
    AND (\"loss of function\" [tiab] OR \"loss-of-function\" [tiab] OR \"haplo-insufficiency\" [tiab] OR haploinsufficiency [tiab] \
    OR \"bi-allelic\" [tiab] OR \"biallelic\" [tiab] OR recessive [tiab] OR homozygous [tiab] OR heterozygous [tiab] OR \"de novo\" \
    [tiab] OR dominant [tiab] OR \" X-linked\" [tiab]) AND (\"intellectual\" [tiab] OR \"mental retardation\" [tiab] OR \"cognitive\" \
    [tiab] OR \"developmental\" [tiab] OR \"neurodevelopmental\" [tiab]) AND “last 2 years”[dp] AND KDM3B)"

    print(Entrez.email)
    query = "((\"2021\"[Date - Publication] : \"3000\"[Date - Publication])) AND (CDH8[Text Word])"


class GeneSearcher:
    def __init__(self):
        self.search_results = None

    def fetch_results(self, parameters: dict) -> int:
        """A method which queries PubMed and returns the count
           of results matching the search parameters. Can be used
           to quickly inform the user if their inserted parameters
           have any results of not.

        :param parameters Parameters entered by the user (dict).
        :return Count of matching articles (int).
        """
        query = parameters.get("input_generated_query")
        min_date: datetime = parameters.get("input_date_after")
        max_date: datetime = parameters.get("input_date_before")
        if query is not None:
            if any(date is not None for date in (min_date, max_date)):
                if min_date is None:
                    raise NoDateAfterSpecified()
                if max_date is None:
                    raise NoDateBeforeSpecified()
                search_results = Entrez.read(
                    Entrez.esearch(
                        db="pubmed",
                        term=query,
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
                        usehistory="y"
                    )
                )
            self.search_results = search_results
            return int(search_results["Count"])
        else:
            raise NoQuerySpecified

    def results_query(self):
        ssl._create_default_https_context = ssl._create_unverified_context
        articles = self.query_pubmed()
        results = []  # results from pubmed en pubtator

        search_results = insert_db()

        for art in articles:
            results.append({"article": {"title": art.article.title,
                                        "pubmed_id": art.article.pubmed_id,
                                        "doi": art.article.doi,
                                        "publication_date": art.article.publication_date,
                                        "abstract": art.article.abstract,
                                        "journal": art.article.journal.name},
                            "genes": art.genes, "diseases": art.diseases})
            search_results.article_list.append({"title": art.article.title,
                                                "pubmed_id": art.article.pubmed_id,
                                                "doi": art.article.doi,
                                                "publication_date": art.article.publication_date,
                                                "abstract": art.article.abstract})
            search_results.journal_list.append({"name": art.article.journal.name})

            for id, gene in art.genes.items():
                if ";" not in id:
                    search_results.genes_list.append(
                        {"ncbi_gene_id": int(id), "hgnc_symbol": str(gene),
                         "in_genepanel": False})
                    search_results.article_gene.append(
                        # Hier moeten de keys de namen van de kolommen bevatten
                        # waar de waardes in moeten komen te staan in de gewenste
                        # tabel.
                        {"article_id": art.article.doi, "gene_id": gene})

        db = DatabaseInserter()
        db.insert_search_results(search_results)
        return results  # list with dict per article

    def query_validator(self, query):
        """
        Validates the query on parentheses being closed off correctly.
        Also checks if search terms contain quotation marks.
        :param query: String which contains the query used as search term
        :return:
        """
        open_list = ["[", "("]
        close_list = ["]", ")"]

        stack = []
        for i in query:
            if i in open_list:
                stack.append(i)
            elif i in close_list:
                pos = close_list.index(i)
                if ((len(stack) > 0) and
                        (open_list[pos] == stack[len(stack) - 1])):
                    stack.pop()
                else:
                    return False
        if len(stack) == 0:
            return True
        else:
            return False

    def query_pubmed(self):
        """
        Uses the query to look for article ids on pubmed.
        Also checks the ids being unique.
        :return:
        """

        count = int(self.search_results["Count"])
        print("Found %i results" % count)

        articles = list()
        batch_size = 100
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

                # print(record)
                pubmed_id = record.get("MedlineCitation").get("PMID")
                title = record.get("MedlineCitation").get("Article").get(
                    "ArticleTitle")

                print(pubmed_id, "pubmed_id")
                doi_element = record.get("MedlineCitation").get("Article").get(
                    "ELocationID")
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
                # print(title, "\n", pubmed_id, "\n", doi, "\n", publication_date,
                #       "\n", abstract, "\n", journal_name)  # example
                art = Article(title=title, pubmed_id=pubmed_id, doi=doi,
                              publication_date=publication_date,
                              abstract=abstract,
                              journal=Journal(name=journal_name))
                # journal = Journal(name=journal_name)
                articles.append(DataArticle(article=art))
                # except (TypeError, KeyError):
                #     print("Element not available")
        # url_maker(articles)
        articles = self.pubtator_output(articles)
        return articles  # articles list with DataArticle complete
        # except None:  # only get None
        #     print(
        #         "The Entrez package is currently offline, please try again later.")

    def extract_date(self, record: dict):
        """A method which safely extracts a publication date
           from an article.

        :param record Article possibly containing a publication date (Dict).
        :return Extracted date from the article if available.
        """

        pub_date = record.get("MedlineCitation").get("Article").get("Journal").get(
            "JournalIssue").get("PubDate")
        if pub_date is not None:
            year = pub_date.get("Year")
            month = pub_date.get("Month")
            day = pub_date.get("Day")
            concatinated = "Not available"
            # x = datetime.datetime(int(year), int(month), int(day))
            if year is not None:
                concatinated = year
                date = datetime.strptime(concatinated, "%Y")
                if month is not None:
                    concatinated += f"-{month}"
                    date = datetime.strptime(concatinated, "%Y-%m")
                if day is not None:
                    concatinated += f"-{day}"
                    print(concatinated)
                    date = datetime.strptime(concatinated, "%Y-%m-%d")
                return date
            return None

    def url_maker(self, idlist):
        """
        Creates the URL used for looking up genes on Pubtator
        :param idlist: List that contains of all the unique ids found in query_pubmed
        :return: complete_url: url for Pubtator
        """

        # post request 1000
        # get request = 100
        url = ""
        for id in idlist:  # id = pubmedid for article pubtator
            if id != idlist[len(idlist) - 1]:
                url += id + ","
            else:
                url += id

        complete_url = "https://www.ncbi.nlm.nih.gov/research/pubtator-api" \
                       "/publications/export/biocxml?pmids={}&" \
                       "concepts=gene,disease".format(url)
        # complete_url is url for pubtator from pubtator API
        return complete_url

    def pubtator_output(self, articles):
        """
        Uses the URL made in url_maker to look up genes on pubtator.
        :param complete_url: string of the url used as input for Pubtator.
        It contains the genes found in the articles
        :return: articles updated with genes from pubtator
        """
        idlist = []
        for article in articles:
            # print(article.article.pubmed_id, "pbid")
            idlist.append(str(article.article.pubmed_id))
        url = self.url_maker(idlist)  # url for all the articles pubmed found
        print(url)
        result = requests.get(url)  # get xml-page pubtator
        if result.status_code == 200:
            print("Request succesful.")
            data = self.parse_results(result)
            if len(data) == len(articles):  # doesn't cont. if something is wrong
                for art in articles:
                    # genes moet waarschijnlijk dict worden ipv genes
                    # print(data)
                    art.genes = data[art.article.pubmed_id][0]  # genes dict
                    art.diseases = data[art.article.pubmed_id][1]  # diseases
                    # print(data[art.article.pubmed_id])

                    # if data[art.article.pubmed_id].keys()[:4] == "MESH":
                    #     art.diseaes
                    # art.genes = data[art.article.pubmed_id]  # added to DataArticle
                    # print(data[art.article.pubmed_id])
                    # not the article object it self yet
            # for a in articles:
            # print(a) # check
        else:
            print("Request not succesful.")
        return articles  # updated DataArticle with genes from pubtator

    def anno_document(self, documents):
        """
        Adds the ncbi gene id and the gene per article to list, so it can
        later be added to the Article object.
        :param documents: <passage>info</passage> all the info from pubtator
        :param count: the amount of articles
        :return:
        """

        data_document = []
        id = 0
        for document in documents.findall("passage"):
            for doc in document.findall("annotation"):
                for anno in doc:
                    if anno.text is None:
                        data_document.append([anno.attrib])

                    elif anno.text == "None":  # when there is no MESH:
                        data_document.append([anno.attrib, "MESH:" + str(id)])
                        id += 1  # if there is no ID
                    else:
                        data_document.append([anno.attrib, anno.text])
        return data_document

    def article_id(self, documents):
        for document in documents.findall("id"):
            return document.text

    def parse_results(self, result):
        """
        Parses all the necessary results out of the Pubtator output.
        :param result: xml form from Pubtator output
        :return: data: list with dict [ncbi gene id, gene]
        """
        tree = etree.fromstring(result.text)
        data_doc = {}
        for documents in tree.findall("document"):
            data_doc[self.article_id(documents)] = self.anno_document(documents)
        data_pubtator = {}
        for pmid, data in data_doc.items():
            # print(pmid, data, "data")
            print(pmid, "pmid")
            genes_pt3 = {}
            mesh = {}
            for gene in data:  # gene = {'key': 'identifier'}, '5362']
                try:
                    if gene[0]["key"] == "identifier":
                        gene_id = gene[1]
                        print(gene_id)
                        genes_pt3[gene_id] = ""
                    if gene[0]["key"] == "Identifier":
                        gene_id = gene[1]
                        mesh[gene_id] = ""
                except KeyError:
                    pass
                if len(gene) >= 2:
                    if not gene[0]:
                        if gene_id[:4] == "MESH":
                            if gene[1].strip() != "":
                                mesh[gene_id.strip()] = gene[1].strip()
                        else:
                            genes_pt3[gene_id] = gene[1]
            data_pubtator[pmid] = [genes_pt3, mesh]
        return data_pubtator

    def article(self, id_list):
        """

        :param id_list: list with pmids f
        :return:
        """

        for article_id in id_list:
            # For finding title, publication date, doi and pubmed ID
            Entrez.email = "Your.Name.Here@example.org"
            handle = Entrez.esummary(db="pubmed", id=article_id)
            record = Entrez.read(handle)
            # print(record)
            article_publication_date = record[0]["EPubDate"]
            article_title = record[0]["Title"]
            article_doi = record[0]["ArticleIds"]["doi"]
            article_pubmed_id = article_id
            handle.close()
            print(article_doi)

            # For finding the abstract
            handle = Entrez.efetch(db="pubmed", id=article_id, rettype="text",
                                   retmode="abstract")
            article_abstract = handle
            handle.close()
            Article1 = Article(title=article_title, pubmed_id=article_pubmed_id,
                               doi=article_doi,
                               publication_date=article_publication_date,
                               abstract=article_abstract)
            print(Article1.title)


# [article, journal, [genes, genes, genes]]

@dataclass
class DataArticle:
    article: Article
    # journal : Journal
    genes: dict = field(default_factory=dict)
    diseases: dict = field(default_factory=dict)


class NoQuerySpecified(Exception):
    """Exception for when no query is specified."""

    def __init__(self):
        super().__init__("No query specified! Can't perform a search when "
                         "no query is specified.")


@dataclass
class insert_db:
    article_list: list = field(default_factory=list)
    genes_list: list = field(default_factory=list)
    journal_list: list = field(default_factory=list)
    article_gene: list = field(default_factory=list)


# def alias_search_hgnc():
#     # # http://rest.genenames.org/fetch/symbol/A2M
#     gene_symbol = "A2M"
#     url = "http://rest.genenames.org/fetch/symbol/{}".format(gene_symbol)
#
#     try:
#         result = requests.get(url)
#         # https://stackoverflow.com/questions/18308529/python-requests-package-handling-xml-response
#         if result.status_code == 200: # webpage legit
#             # handle = Entrez.parse(result.text)
#             # record = Entrez.read(handle)
#             # print(record)
#             tree = ElementTree.fromstring(result.text)
#             print(tree)
#             for child in tree:
#                 print(child)
#             # root = ET.fromstring(country_data_as_string)
#
#             # parse hier de xml output
#
#     except requests.exceptions.RequestException:
#         print("Geen entry gevonden op HGNC")

class NoDateAfterSpecified(Exception):
    """Exception for when parameter 'date after' is not specified
       but 'date before' is specified.
    """

    def __init__(self):
        super().__init__("No 'date after' specified! 'Date after'-parameter "
                         "is required when 'date before'-parameter is specified.")


class NoDateBeforeSpecified(Exception):
    """Exception for when parameter 'date before' is not specified
       but 'date after' is specified.
    """

    def __init__(self):
        super().__init__("No 'date before' specified! 'Date before'-parameter "
                         "is required when 'date after'-parameter is specified.")


if __name__ == "__main__":
    main()
