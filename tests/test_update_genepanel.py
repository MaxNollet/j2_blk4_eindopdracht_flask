import ntpath
from os import path
from typing import Union

from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support import expected_conditions
from selenium.webdriver.support.ui import WebDriverWait

from tests.test_query_builder import HelperFunctions

webDriver = Union[
    webdriver.Firefox,
    webdriver.Chrome,
    webdriver.Edge,
    webdriver.Safari,
    webdriver.Ie,
    webdriver.Opera,
    webdriver.PhantomJS
]


class ElementSelectorUpdateGenepanel:
    """a class that groups methods which select various
       elements on the update_genepanel-page.

    base_url = URL of the page containing the elements (str).
    """
    base_url = "http://127.0.0.1:5000/update_genepanel"

    @staticmethod
    def select_input_upload_genepanel(selenium: webDriver):
        """Slect element 'input_upload_genepanel'."""
        return selenium.find_element_by_id("input_upload_genepanel")

    @staticmethod
    def select_button_clear(selenium: webDriver):
        """Select element 'button_clear'."""
        return selenium.find_element_by_id("button_clear")

    @staticmethod
    def select_button_update_genepanel(selenium: webDriver):
        """Select element 'button_update_genepanel'."""
        return selenium.find_element_by_id("button_update_genepanel")


class TestDefaultsUpdateGenepanel(ElementSelectorUpdateGenepanel):
    """A class which groups tests that check the defaults
       for the update_genepanel-page.
    """

    def test_input_upload_genepanel_enabled(self, selenium: webDriver):
        """Test input_upload_genepanel if enabled."""
        selenium.get(self.base_url)
        input_upload_genepanel = self.select_input_upload_genepanel(selenium)
        assert input_upload_genepanel.is_enabled() is True

    def test_input_upload_genepanel_default_value(self, selenium: webDriver):
        """Test input_upload_genepanel default value."""
        selenium.get(self.base_url)
        input_upload_genepanel = self.select_input_upload_genepanel(selenium)
        default_value = input_upload_genepanel.get_attribute("value")
        assert default_value == ""

    def test_button_clear_enabled(self, selenium: webDriver):
        """Test button_clear if enabled."""
        selenium.get(self.base_url)
        button_clear = self.select_button_clear(selenium)
        assert button_clear.is_enabled() is True

    def test_button_clear_text(self, selenium: webDriver):
        """Test button_clear if text is correct."""
        selenium.get(self.base_url)
        button_clear = self.select_button_clear(selenium)
        assert button_clear.text == "Clear"

    def test_button_update_genepanel_enabled(self, selenium: webDriver):
        """Test button_update_genepanel is enabled."""
        selenium.get(self.base_url)
        button_update_genepanel = self.select_button_update_genepanel(selenium)
        assert button_update_genepanel.is_enabled() is True

    def test_button_update_genepanel_type(self, selenium: webDriver):
        """Test button_update_genepanel if correct type."""
        selenium.get(self.base_url)
        button_update_genepanel = self.select_button_update_genepanel(selenium)
        button_type = button_update_genepanel.get_attribute("type")
        assert button_type == "submit"


# class TestJavaScript(ElementSelectorUpdateGenepanel, HelperFunctions):
#     """A class that groups tests which are related to
#        JavaScript-based features.
#     """
#
#     def test_select_file(self, selenium: webDriver):
#         """Test selecting a file."""
#         selenium.get(self.base_url)
#         input_upload_genepanel = self.select_input_upload_genepanel(selenium)
#



# class TestUpdateGenepanel(ElementSelectorUpdateGenepanel):
#     """A class which groups tests that check if a genepanel
#        can be updated.
#     """
#
#     def test_upload_file(self, selenium: webDriver):
#         """Test for uploading files"""
#         selenium.get(self.base_url)
#



# def test_clear_file(selenium: webDriver):
#     """Test the button that clears a file form the file-chooser."""
#     selenium.get("http://127.0.0.1:5000/update_genepanel")
#     input_upload_genepanel = selenium.find_element_by_id("input_upload_genepanel")
#     button_clear = selenium.find_element_by_id("button_clear")
#     file_path = path.abspath(path.join(path.dirname(__file__), "..", "requirements.txt"))
#     input_upload_genepanel.send_keys(file_path)
#     assert ntpath.basename(input_upload_genepanel.get_attribute("value")) == ntpath.basename(file_path)
#     WebDriverWait(selenium, 10).until(
#         expected_conditions.element_to_be_clickable((By.ID, "input_upload_genepanel"))
#     )
#     selenium.execute_script("arguments[0].click();", button_clear)
#     assert input_upload_genepanel.get_attribute("value") == ""
