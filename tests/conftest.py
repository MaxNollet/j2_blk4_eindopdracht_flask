import os
import subprocess

import pytest


def pytest_addoption(parser):
    """A function which adds a new option to PyTest:
       the --headless option makes sure to run Selenium-
       tests in a headless-mode when the browser supports
       the option.
    """
    group = parser.getgroup('selenium', 'selenium')
    group.addoption('--headless',
                    action='store_true',
                    help='enable headless mode for supported browsers.')


@pytest.fixture
def chrome_options(chrome_options, pytestconfig):
    """A function which passes the 'headless'-option
       to Chrome to run Selenium-tests in a headless-
       mode.
    """
    if pytestconfig.getoption('headless'):
        chrome_options.add_argument('headless')
    return chrome_options


@pytest.fixture(scope="session", autouse=True)
def live_server():
    """A function which tries to solve an issue when the
       original liver_server-fixture tries to pickle an
       object but fails on Windows 10.
    """
    env = os.environ.copy()
    env["FLASK_APP"] = "gaps"
    server = subprocess.Popen(['flask', 'run', '--port', "5000"], env=env)
    try:
        yield server
    finally:
        server.terminate()
