from flask import Blueprint, render_template, request, url_for

home_page = Blueprint("home_page", __name__)


@home_page.route("/")
def homepage():
    return render_template("homepage.html", active="home")

# shift command r = hard refresh cache/ en zooi
