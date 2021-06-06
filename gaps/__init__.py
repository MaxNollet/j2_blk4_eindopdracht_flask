import os
from os import environ

from flask import Flask

from gaps.blueprint_api import blueprint_api
from gaps.blueprint_homepage import blueprint_homepage
from gaps.blueprint_query_builder import blueprint_query_builder
from gaps.blueprint_results import blueprint_results
from gaps.blueprint_update_genepanel import blueprint_update_genepanel
from gaps.models import db


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
            app.config.from_object("config.Production")
        else:
            value = environ.get('FLASK_ENV').lower().capitalize()
            try:
                app.config.from_object(f"config.{value}")
            except ImportError:
                app.config.from_object("config.Production")
    else:
        app.config.from_object("config.Testing")

    # Prepare the upload-folder for file uploads.
    upload_path = app.config['UPLOAD_PATH']
    if not os.path.exists(upload_path):
        os.mkdir(upload_path)
    else:
        for file in os.listdir(upload_path):
            os.remove(os.path.join(upload_path, file))

    # Initiate database and register blueprints.
    db.init_app(app)
    app.register_blueprint(blueprint_homepage)
    app.register_blueprint(blueprint_query_builder)
    app.register_blueprint(blueprint_update_genepanel)
    app.register_blueprint(blueprint_api)
    app.register_blueprint(blueprint_results)
    return app
