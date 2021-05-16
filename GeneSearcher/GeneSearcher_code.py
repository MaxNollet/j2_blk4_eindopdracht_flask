from Bio import Entrez
import requests
import pyhgnc

def main():
    query = "((variant [tiab] OR variants [tiab] OR mutation [tiab] OR mutations [tiab] OR substitutions [tiab] OR substitution [tiab] ) \
    AND (\"loss of function\" [tiab] OR \"loss-of-function\" [tiab] OR \"haplo-insufficiency\" [tiab] OR haploinsufficiency [tiab] \
    OR \"bi-allelic\" [tiab] OR \"biallelic\" [tiab] OR recessive [tiab] OR homozygous [tiab] OR heterozygous [tiab] OR \"de novo\" \
    [tiab] OR dominant [tiab] OR \" X-linked\" [tiab]) AND (\"intellectual\" [tiab] OR \"mental retardation\" [tiab] OR \"cognitive\" \
    [tiab] OR \"developmental\" [tiab] OR \"neurodevelopmental\" [tiab]) AND “last 2 years”[dp] AND KDM3B)"

    query2 = "((\"2021\"[Date - Publication] : \"3000\"[Date - Publication])) AND (CDH8[Text Word])"

    idlist = query_pubmed(query)
    hele_url = url_maker(idlist)
    pubtator_output(hele_url)
    #query_HGNC("AGPAT4")

def query_validator(query):
    """
    Validates the query on parentheses being closed off correctly.
    Also checks if search terms contain quotation marks.
    :param query: String which contains the query used as search term
    :return: 
    """
    open_list = ["[","("]
    close_list = ["]",")"]

    stack = []
    for i in query:
        if i in open_list:
            stack.append(i)
        elif i in close_list:
            pos = close_list.index(i)
            if ((len(stack) > 0) and
                    (open_list[pos] == stack[len(stack)-1])):
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
    query_validator(query)
    if query_validator(query):
        print("The query is valid.")
        Entrez.email = "femke.nijman@outlook.com"
        searchhandle = Entrez.esearch(db="pubmed", term=query)
        search_results = Entrez.read(searchhandle)

        idlist = search_results["IdList"]
        if(len(set(idlist)) == len(idlist)):
            list_unique = True
            print("List does not contain duplicates.")
        else:
            list_unique = False
            print("List does not contain duplicates")
        print(idlist)
        return idlist
    else:
        print("The query isn't valid")

def url_maker(idlist):
    """
    Creates the URL used for looking up genes on Pubtator
    :param idlist: List that contains of all the unique ids found in query_pubmed
    :return: 
    """
    url = ""
    for i in idlist:
        if i != idlist[len(idlist)-1]:
            url += i + ","
        else:
            url += i
    complete_url = "https://www.ncbi.nlm.nih.gov/research/pubtator-api/publications/export/pubtator?pmids=" + url + "&concepts=gene"
    print(complete_url)
    return complete_url

def pubtator_output(complete_url):
    """
    Uses the URL made in url_maker to look up genes on pubtator.
    :param complete_url: string of the url used as input for Pubtator. It contains the genes found in the articles
    :return:
    """
    result = requests.get(complete_url)
    status_code = result.status_code
    if status_code == 200:
        print("Request succesful.")
    else:
        print("Request not succesful.")
    #print(result.text)

#def query_HGNC(gene):
    #pyhgnc.set_mysql_connection(host='localhost', user='pyhgnc_user', passwd='pyhgnc_passwd', db='pyhgnc')
    #query = pyhgnc.query()
    #result = query.alias_symbol(alias_symbol=gene)[0]
    #print(result.hgnc)

main()
