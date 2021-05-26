from flask import Blueprint, request, jsonify, render_template
from typing import Tuple, Dict
from werkzeug.utils import secure_filename
from gaps.gene_searcher import GeneSearcher_code as retriever

blueprint_api = Blueprint("blueprint_api", __name__)


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
    # main()
    secure_filenames = VerifyFormParameters.get_valid_filenames(
        "input_load_symbols")
    valid_parameters = VerifyFormParameters.get_valid_parameters()
    # genesearch_code(valid_parameters["input_generated_query"])
    response = {
        "input": {"files": secure_filenames, "parameters": valid_parameters}}
    results = retriever.results_query(
        valid_parameters["input_generated_query"])
    print(results)

    return render_template("template_results.html", results=results)
    # return jsonify(response)


@blueprint_api.route("/update_genepanel", methods=["POST"])
def update_genepanel_submit():
    """A function which handles POST-requests to the
       '/update_genepanel'-route. This route is used
       to inform the client about valid uploaded
       files and what has changed after the update of
       a genepanel.

    :return JSON-response containing valid values (JSON).
    """
    secure_filenames = VerifyFormParameters.get_valid_filenames("input_upload_genepanel")
    response = {"input": {"files": secure_filenames}}
    return jsonify(response)


class VerifyFormParameters:
    """A class that extracts valid/secure variables and
       from requests so the server can safely process
       these variables and files.
    """

    @staticmethod
    def get_valid_filenames(file_upload_field: str) -> Tuple[str]:
        """A method which secures filenames from requests
           and returns the secured filenames for safe
           storage and processing of the files.

        :param file_upload_field Field to check containing
                                 possible filenames (str).
        :return Secured filenames (tuple(str)).
        """
        if file_upload_field in request.files:
            uploaded_files = request.files.getlist(file_upload_field)
            filenames = list()
            for uploaded_file in uploaded_files:
                filename = secure_filename(uploaded_file.filename)
                if filename != "":
                    filenames.append(filename)
            if len(filenames) > 0:
                return tuple(filenames)

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
                verified_values[parameter_key] = submitted_value
        if len(verified_values) > 0:
            return verified_values
