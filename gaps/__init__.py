from os import environ

from flask import Flask

from gaps.home_page import home_page
from gaps.query import query_page
from gaps.models import db


def create_app():
    """A function which configures the application and
       the database. Uses the environment-variable
       'FLASK_ENV' to choose a config-profile.

    :return Configured Flask-app (Flask).
    """
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
    return app
