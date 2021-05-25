class Article:
    def __init__(self, title, pubmed_id, doi, publication_date, abstract):
        self.title = title
        self.pubmed_id = pubmed_id
        self.doi = doi
        self.publication_date = publication_date
        self.abstract = abstract
      #  self.journal_id = journal_id

    def get_title(self):
        return self.title

    def get_pubmed_id(self):
        return self.pubmed_id

    def get_doi(self):
        return self.doi

    def get_publication_date(self):
        return self.publication_date

    def get_abstract(self):
        return self.abstract
