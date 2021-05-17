import pytest
import os
import sys

current_dir = os.path.abspath(os.path.dirname(__file__))
parent_dir = os.path.dirname(current_dir)
sys.path.append(parent_dir)

from gaps import create_app
from gaps.models import db


@pytest.fixture
def client():
    """A function which prepares the application for
       testing and returns a testing client.

    :return Testing client (Flask).
    """
    app = create_app(testing=True)
    with app.app_context():
        db.create_all()
    yield app.test_client()
