from flask import Flask, render_template
from .home_page import home_page    # . voor de structuur
from .query import query_page


def create_app(test_config=None):
    app = Flask(__name__, instance_relative_config=True)
    app.config.from_mapping(SECRET_KEY="dev")

    if test_config is None:
        app.config.from_pyfile("config.py", silent=True)
    else:
        app.config.from_mapping(test_config)

    app.register_blueprint(home_page)

    app.register_blueprint(query_page)

    # @app.route("/")
    # def hello_world():
    #     return render_template("hello_flask.html")

    return app
