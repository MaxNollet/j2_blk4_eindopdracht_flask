import os.path
import random
from datetime import datetime
from typing import Dict

from flask import Blueprint, request, jsonify, current_app, abort, url_for
from werkzeug.datastructures import FileStorage
from werkzeug.utils import secure_filename

from gaps.gene_searcher import *
from gaps import file_reader as reader
from gaps.genelogic import GenepanelReader, DatabaseInserter, \
    GenepanelColumnNotFound

blueprint_api = Blueprint("blueprint_api", __name__)


def genes_box(searcher, valid_parameters):
    if valid_parameters["radio_include_symbols"] == "true":
        searcher.include_exclude = True
    else:
        searcher.include_exclude = False
    genes_list = list(valid_parameters["input_symbols"].split(","))
    set_genes_list = set()
    for gene in genes_list:
        set_genes_list.add(gene.strip())
    searcher.specific_gene_symbols = set_genes_list


@blueprint_api.route("/query_builder", methods=["POST"])
def query_builder_submit():
    """A function which handles POST-requests to the
       '/query_builder'-route. This route is used to
       inform the client about uploaded valid files,
       if the query has any results and to redirect
       the user to the results-page if there are
       results.

    :return JSON-response containing valid values (JSON).
    """
    secure_filenames = VerifyFormParameters.get_valid_filenames(
        "input_load_symbols")
    valid_parameters = VerifyFormParameters.get_valid_parameters()
    upload_path = current_app.config['UPLOAD_PATH']
    try:
        if secure_filenames:
            filenames = tuple(secure_filenames.keys())
            filename = filenames[0]
            file_location = os.path.join(upload_path, filename)
            secure_filenames[filename].save(file_location)
        print(valid_parameters)
        print(secure_filenames)
        searcher = GeneSearcher()
        # check dit met true of false

        if valid_parameters.get("input_symbols"):
            searcher.specify_gene = True
            genes_box(searcher, valid_parameters)
        elif file_location and len(file_location) > 1:  # for the file input
            searcher.specify_gene = True  # use white/black list
            if valid_parameters["radio_include_symbols"] == "true":
                searcher.include_exclude = True
            else:
                searcher.include_exclude = False
            searcher.specific_gene_symbols = reader.file_reader(file_location)

            print("Reading file")
        elif file_location and len(file_location) > 1 and \
                valid_parameters.get("input_symbols"):
            print("komt die")
            searcher.specify_gene = True
            genes_box(searcher, valid_parameters)
        count = searcher.fetch_results(valid_parameters)
        if count < 1:
            response = {
                "message": "No articles found! Adjust your search parameters "
                           "and try again.",
                "type": "info"}
        else:
            searcher.results_query()
            response = {"redirect": url_for("blueprint_results.results",
                                            query_id=searcher.uuid_query)}
        return jsonify(response)
    except (NoQuerySpecified, NoDateAfterSpecified, NoDateBeforeSpecified,
            MalformedQuery, IncorrectArticleFound, NoGeneFound) as e:
        response = {"message": str(e),
                    "type": "warning"}
    finally:
        if secure_filenames:
            for filename in secure_filenames.keys():
                os.remove(os.path.join(upload_path, filename))
    return jsonify(response)


@blueprint_api.route("/update_genepanel", methods=["POST"])
def update_genepanel_submit():
    """A function which handles POST-requests to the
       '/update_genepanel'-route. This route is used
       to inform the client about valid uploaded
       files and what has changed after the update of
       a genepanel.

    :return JSON-response containing valid values (JSON).
    """
    secure_filenames = VerifyFormParameters.get_valid_filenames(
        "input_upload_genepanel")
    upload_path = current_app.config['UPLOAD_PATH']
    try:
        if secure_filenames:
            filenames = tuple(secure_filenames.keys())
            filename = filenames[0]
            file_location = os.path.join(upload_path, filename)
            secure_filenames[filename].save(file_location)
            genepanel_data = GenepanelReader(
                filename=file_location).get_genepanel_data()
            DatabaseInserter().insert_genepanel(genepanel_data)
            response = {
                "message": "Genepanels updated successfully! Refresh the "
                           "page to see the new statistics about "
                           "the updated geneanels.", "type": "success"}
        else:
            response = {
                "message": "No file uploaded! Could not update genepanels.",
                "type": "danger"}
    except GenepanelColumnNotFound as e:
        response = {"message": str(e),
                    "type": "danger"}
    finally:
        if secure_filenames:
            for filename in secure_filenames.keys():
                os.remove(os.path.join(upload_path, filename))
    return jsonify(response)


class VerifyFormParameters:
    """A class that extracts valid/secure variables and
       from requests so the server can safely process
       these variables and files.
    """

    @staticmethod
    def get_valid_filenames(file_upload_field: str) -> Dict[str, FileStorage]:
        """A method which secures filenames from requests
           and returns the secured filenames for safe
           storage and processing of the files.

        :param file_upload_field Field to check containing
                                 possible filenames (str).
        :return Secured filenames (Dict[str, FileStorage]).
        """
        if file_upload_field in request.files:
            uploaded_files = request.files.getlist(file_upload_field)
            filenames = dict()
            for uploaded_file in uploaded_files:
                if uploaded_file.filename != "":
                    file_extension = os.path.splitext(uploaded_file.filename)[
                        1]
                    if file_extension in current_app.config[
                        "UPLOAD_EXTENSIONS"]:
                        filename = secure_filename(uploaded_file.filename)
                        if filename != "":
                            filenames[
                                f"{random.randint(0, 9000)}_{filename}"] = uploaded_file
                    else:
                        abort(400)
            if len(filenames) > 0:
                return filenames

    @staticmethod
    def get_valid_parameters() -> Dict[str, str]:
        """A method which extracts all elements from a
           form that are not empty and returns them for
           further processing.

        :return All non-default values from the form
                (dict(str, str)).
        """
        verified_values = dict()
        for parameter_key in request.form.keys():
            submitted_value = request.form[parameter_key].strip()
            if submitted_value != "":
                if "date" in parameter_key:
                    submitted_value = datetime.strptime(submitted_value,
                                                        "%Y-%m-%d")
                verified_values[parameter_key] = submitted_value
        if len(verified_values) > 0:
            return verified_values
