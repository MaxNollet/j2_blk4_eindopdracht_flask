import ntpath
from os import path
from typing import Union

from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support import expected_conditions
from selenium.webdriver.support.ui import WebDriverWait

webDriver = Union[
    webdriver.Firefox,
    webdriver.Chrome,
    webdriver.Edge,
    webdriver.Safari,
    webdriver.Ie,
    webdriver.Opera,
    webdriver.PhantomJS
]


def test_clear_file(selenium: webDriver):
    """Test the button that clears a file form the file-chooser."""
    selenium.get("http://127.0.0.1:5000/update_genepanel")
    input_upload_genepanel = selenium.find_element_by_id("input_upload_genepanel")
    button_clear = selenium.find_element_by_id("button_clear")
    file_path = path.abspath(path.join(path.dirname(__file__), "..", "requirements.txt"))
    input_upload_genepanel.send_keys(file_path)
    assert ntpath.basename(input_upload_genepanel.get_attribute("value")) == ntpath.basename(file_path)
    WebDriverWait(selenium, 10).until(
        expected_conditions.element_to_be_clickable((By.ID, "input_upload_genepanel"))
    )
    selenium.execute_script("arguments[0].click();", button_clear)
    assert input_upload_genepanel.get_attribute("value") == ""
