from flask import Blueprint, render_template

query_page = Blueprint("query_page", __name__)


@query_page.route("/query")
def homepage():
    return render_template("query.html", active="query_input")

# shift command r = hard refresh cache/ en zooi