from os import environ, path
from dotenv import load_dotenv

base_directory = path.abspath(path.dirname(__file__))
load_dotenv(path.join(base_directory, ".env"))


class Config(object):
    SECRET_KEY = environ.get("SECRET_KEY")
    # Config for the database
    SQLALCHEMY_DATABASE_URI = environ.get("SQLALCHEMY_DATABASE_URI")
    SQLALCHEMY_TRACK_MODIFICATIONS = False


class Production(Config):
    FLASK_ENV = "production"
    DEBUG = False
    TESTING = False
    # Config for the database
    SQLALCHEMY_ECHO = False


class Development(Config):
    FLASK_ENV = "development"
    DEBUG = True
    TESTING = True
    # Config for the database
    SQLALCHEMY_ECHO = True
