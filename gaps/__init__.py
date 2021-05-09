from os import environ

from flask import Flask

from gaps.home_page import home_page  # . voor de structuur
from gaps.query import query_page
from gaps.models import db


def create_app():
    app = Flask(__name__, instance_relative_config=True)
    if environ.get("FLASK_ENV") is None:
        app.config.from_object("config.Production", silent=True)
    else:
        value = environ.get('FLASK_ENV').lower().capitalize()
        try:
            app.config.from_object(f"config.{value}")
        except ImportError:
            app.config.from_object("config.Production", silent=True)
    db.init_app(app)
    app.register_blueprint(home_page)
    app.register_blueprint(query_page)

    # @app.route("/")
    # def hello_world():
    #     return render_template("hello_flask.html")

    return app
