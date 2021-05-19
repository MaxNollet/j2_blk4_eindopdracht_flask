from os import environ

from flask import Flask

from gaps.blueprint_home import blueprint_home
from gaps.blueprint_query import blueprint_query
from gaps.models import db, Gene
from gaps.genelogic import reader, DatabaseInserter


def create_app(testing=False):
    """A function which configures the application and
       the database. Uses the environment-variable
       'FLASK_ENV' to choose a config-profile.

    :param testing Load the testing-configuration or not (bool).
    :return Configured Flask-app (Flask).
    """
    app = Flask(__name__, instance_relative_config=True)
    if not testing:
        if environ.get("FLASK_ENV") is None:
            app.config.from_object("config.Production", silent=True)
        else:
            value = environ.get('FLASK_ENV').lower().capitalize()
            try:
                app.config.from_object(f"config.{value}")
            except ImportError:
                app.config.from_object("config.Production")
    else:
        app.config.from_object("config.Testing")
    # Initiate database and register blueprints.
    DatabaseInserter.updateGenpanel(db)
    # print(db)
    db.init_app(app)

    app.register_blueprint(blueprint_home)
    app.register_blueprint(blueprint_query)
    return app
