from dataclasses import dataclass, field
from xml.etree import ElementTree as etree
from os import environ
import requests
import ssl
from Bio import Entrez
from gaps.models import Article, Journal, Gene

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

    articles = query_pubmed(query)
    # print(articles)
    # hele_url = url_maker(["33833667", "33810959"])
    # pubtator_output(hele_url)
    # query_HGNC("AGPAT4")
    # alias_search_hgnc()
    # article(idlist)


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


def query_pubmed(query):
    """
    Uses the query to look for article ids on pubmed.
    Also checks the ids being unique. 
    :param query: String which contains the query used as search term
    :return:
    """
    search_results = Entrez.read(
        Entrez.esearch(
            db="pubmed",
            term=query,
            usehistory="y"
        )
    )
    count = int(search_results["Count"])
    print("Found %i results" % count)

    articles = list()
    batch_size = 10
    for start in range(0, count, batch_size):
        handle = Entrez.efetch(db="pubmed",
                               rettype="medline",
                               retmode="xml",
                               retstart=start,
                               retmax=batch_size,
                               webenv=search_results["WebEnv"],
                               query_key=search_results["QueryKey"],
                               )
        records = Entrez.read(handle)
        for record in records["PubmedArticle"]:
            # print(record)
            title = record["MedlineCitation"]["Article"]["ArticleTitle"]
            pubmed_id = record["MedlineCitation"]["PMID"]
            print(pubmed_id, "pubmed_id")
            doi = record["MedlineCitation"]["Article"]["ELocationID"][0]
            # publication_year = record["MedlineCitation"]["Article"]["ArticleDate"]["Year"]
            publication_year = \
                record["MedlineCitation"]["Article"]["ArticleDate"][0]["Year"]
            publication_month = \
                record["MedlineCitation"]["Article"]["ArticleDate"][0]["Month"]
            publication_day = \
                record["MedlineCitation"]["Article"]["ArticleDate"][0]["Day"]
            # jaar maand dag
            publication_date = publication_year + "-" + \
                               publication_month + "-" + publication_day
            abstract = \
                record["MedlineCitation"]["Article"]["Abstract"][
                    "AbstractText"][0]
            journal_name = record["MedlineCitation"]["Article"]["Journal"][
                "Title"]
            # print(title, "\n", pubmed_id, "\n", doi, "\n", publication_date,
            #       "\n", abstract, "\n", journal_name)  # example
            art = Article(title=title, pubmed_id=pubmed_id, doi=doi,
                          publication_date=publication_date, abstract=abstract,
                          journal=Journal(name=journal_name))
            # journal = Journal(name=journal_name)
            articles.append(DataArticle(article=art))
    # url_maker(articles)
    articles = pubtator_output(articles)
    return articles  # articles list with DataArticle complete

def url_maker(idlist):
    """
    Creates the URL used for looking up genes on Pubtator
    :param idlist: List that contains of all the unique ids found in query_pubmed
    :return: complete_url: url for Pubtator
    """

    # print(article.article.pubmed_id, "pbid")
    url = ""
    for id in idlist:  # id = pubmedid for article pubtator
        if id != idlist[len(idlist) - 1]:
            url += id + ","
        else:
            url += id

    complete_url = "https://www.ncbi.nlm.nih.gov/research/pubtator-api" \
                   "/publications/export/biocxml?pmids={}&" \
                   "concepts=gene".format(url)
    # complete_url is url for pubtator from pubtator API
    return complete_url


def pubtator_output(articles):
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
    url = url_maker(idlist)  # url for all the articles pubmed found
    result = requests.get(url)  # get xml-page pubtator
    if result.status_code == 200:
        print("Request succesful.")
        data = parse_results(result)
        if len(data) == len(articles):  # so they can use the same index
            for index, article in enumerate(articles):
                # genes moet waarschijnlijk dict worden ipv genes
                print(index)
                article.genes = data[index]  # added to DataArticle
                # not the article object it self yet

        # for a in articles:
        # print(a) # check

    else:
        print("Request not succesful.")
    return articles  # updated DataArticle with genes from pubtator


def anno_document(documents, count):
    """
    Adds the ncbi gene id and the gene per article to list, so it can
    later be added to the Article object.
    :param documents: <passage>info</passage> all the info from pubtator
    :param count: the amount of articles
    :return:
    """
    data_document = []
    for document in documents.findall("passage"):
        for doc in document.findall("annotation"):
            for anno in doc:
                data_document.append([count, anno.attrib, anno.text])
    return data_document


def parse_results(result):
    """
    Parses all the necessary results out of the Pubtator output.
    :param result: xml form from Pubtator output
    :return: data: list with dict [ncbi gene id, gene]
    """
    tree = etree.fromstring(result.text)
    data_documents = []
    count = -1
    for documents in tree.findall("document"):
        count += 1
        data_documents.append(anno_document(documents, count))

    data = []
    data2 = []  # contains only the unique ncbi gene id / genes from pubtator
    c = 0
    for gene_idlist in data_documents:
        # print(t, "hats")
        genes_pt = {}  # contains all the genes in the article from pubtator
        genes_pt2 = {}  # contains only the unique ncbi gene id / genes from pubtator

        for gi in gene_idlist:
            try:
                if gi[1]["key"] == "identifier":
                    c += 1  # as key because the ncbi_gene id is not unique in article
                    genes_pt[c] = gi[2]  # ncbi gene id
                    iden = gi[2]
                    genes_pt2[iden] = ""
            except KeyError:
                pass  # only need the identifier key
            if not gi[1]:  # de annotion dict/ text is always empty
                genes_pt[c] = [genes_pt[c], gi[2]]
                genes_pt2[iden] = gi[2]
        data.append(genes_pt)
        data2.append(genes_pt2)
    print(data2)
    print(len(data[0]))
    print(len(data2[0]))
    print(len(data[1]))
    print(len(data2[1]))
    return data2  # onlt the unique genes


def article(id_list):
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


# http://rest.genenames.org/fetch/symbol/A2M
# gene_symbol # wordt het symbol voor in de link


# def query_HGNC(gene):
# pyhgnc.set_mysql_connection(host='localhost', user='pyhgnc_user', passwd='pyhgnc_passwd', db='pyhgnc')
# query = pyhgnc.query()
# result = query.alias_symbol(alias_symbol=gene)[0]
# print(result.hgnc)


if __name__ == "__main__":
    main()
