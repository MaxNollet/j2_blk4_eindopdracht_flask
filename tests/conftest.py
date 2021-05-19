import os
import socket
import subprocess

import pytest


@pytest.fixture(scope="session")
def flask_port():
    # Ask OS for a free port.
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.bind(("", 0))
        addr = s.getsockname()
        port = addr[1]
        return port


@pytest.fixture(scope="session", autouse=True)
def live_server(flask_port):
    env = os.environ.copy()
    env["FLASK_APP"] = "gaps"
    # server = subprocess.Popen(['flask', 'run', '--port', str(flask_port)], env=env)
    server = subprocess.Popen(['flask', 'run', '--port', "5000"], env=env)
    try:
        yield server
    finally:
        server.terminate()
