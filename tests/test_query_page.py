from typing import Union

from selenium import webdriver
from selenium.webdriver.support.select import Select

webDriver = Union[
    webdriver.Firefox,
    webdriver.Chrome,
    webdriver.Edge,
    webdriver.Safari,
    webdriver.Ie,
    webdriver.Opera,
    webdriver.PhantomJS
]


class TestDefaults:
    """A class which groups tests related to checking
       default values and settings for elements in the
       interface.
    """

    def test_defaults_query_generator(self, selenium: webDriver, base_url):
        """Test the default values for the query generator."""
        selenium.get("http://127.0.0.1:5000/")
        input_field = selenium.find_element_by_id("input_field")
        input_search_term = selenium.find_element_by_id("input_search_term")
        input_add_type = selenium.find_element_by_id("input_add_type")
        input_generated_query = selenium.find_element_by_id("input_generated_query")
        select_input_field = Select(input_field)
        select_input_add_type = Select(input_add_type)
        # Assert defaults for input_field.
        assert input_field.is_enabled() is True
        assert select_input_field.first_selected_option.text == "All Fields"
        assert len(select_input_field.options) == 51
        # Asserts defaults for input_search_term.
        assert input_search_term.is_enabled() is True
        assert input_search_term.get_attribute("value") == ""
        # Assert defaults for input_add_type.
        assert input_add_type.is_enabled() is False
        assert select_input_add_type.first_selected_option.text == "AND"
        assert len(select_input_add_type.options) == 3
        # Assert defaults for input_generated_query.
        assert input_generated_query.is_enabled() is True
        assert input_generated_query.get_attribute("value") == ""

    def test_defaults_gene_symbols(self, selenium: webDriver):
        """Test the default values for specifying genes."""
        selenium.get("http://127.0.0.1:5000/")
        radio_include_symbols = selenium.find_element_by_id("radio_include_symbols")
        radio_exclude_symbols = selenium.find_element_by_id("radio_exclude_symbols")
        input_symbols = selenium.find_element_by_id("input_symbols")
        input_load_symbols = selenium.find_element_by_id("input_load_symbols")
        button_clear_file = selenium.find_element_by_id("button_clear_file")
        # Assert defaults for radio-buttons.
        assert radio_include_symbols.is_enabled() is True
        assert radio_include_symbols.is_selected() is True
        assert radio_include_symbols.get_attribute("value") == "true"
        assert radio_exclude_symbols.is_enabled() is True
        assert radio_exclude_symbols.is_selected() is False
        assert radio_exclude_symbols.get_attribute("value") == "false"
        # Assert defaults for gene symbol input.
        assert input_symbols.is_enabled() is True
        assert input_symbols.get_attribute("value") == ""
        # Assert defaults for file chooser.
        assert input_load_symbols.is_enabled() is True
        assert input_load_symbols.get_attribute("value") == ""
        # Assert defaults for clear-button.
        assert button_clear_file.is_enabled() is True

    def test_defaults_optional_options(self, selenium: webDriver):
        """Test the default values for optional options and submit-buttons."""
        selenium.get("http://127.0.0.1:5000/")
        input_date_after = selenium.find_element_by_id("input_date_after")
        input_date_before = selenium.find_element_by_id("input_date_before")
        check_new_tab = selenium.find_element_by_id("check_new_tab")
        button_submit = selenium.find_element_by_id("button_submit")
        button_clear = selenium.find_element_by_id("button_clear")
        # Assert defaults for date-inputs.
        assert input_date_after.is_enabled() is True
        assert input_date_after.get_attribute("value") == ""
        assert input_date_before.is_enabled() is True
        assert input_date_before.get_attribute("value") == ""
        # Assert defaults for check open in new tab, submit- and rest-buttons.
        assert check_new_tab.is_enabled() is True
        assert check_new_tab.is_selected() is False
        assert button_submit.is_enabled() is True
        assert button_submit.get_attribute("type") == "submit"
        assert button_clear.is_enabled() is True
        assert button_clear.get_attribute("type") == "reset"


class TestQueryBuilder:
    """A class which groups tests related to the query
       builder. These tests are in place to ensure the
       query builder keeps working as intended.
    """

    def test_first_addition_no_spaces_button(self, selenium: webDriver):
        """Test the addition of a term to the query without spaces."""
        selenium.get("http://127.0.0.1:5000/")
        input_term = selenium.find_element_by_id("input_search_term")
        input_query = selenium.find_element_by_id("input_generated_query")
        button_add = selenium.find_element_by_id("button_add_item")
        input_term.send_keys("Alzheimer")
        button_add.click()
        assert input_query.get_attribute("value") == "Alzheimer[ALL]"

    def test_first_addition_with_spaces_button(self, selenium: webDriver):
        """Test the addition of a term to the query with spaces."""
        selenium.get("http://127.0.0.1:5000/")
        input_term = selenium.find_element_by_id("input_search_term")
        input_query = selenium.find_element_by_id("input_generated_query")
        button_add = selenium.find_element_by_id("button_add_item")
        input_term.send_keys("Frikandel speciaal")
        button_add.click()
        assert input_query.get_attribute("value") == "Frikandel speciaal[ALL]"


class TestJavaScript:
    """A class which groups tests related to JavaScript-
       functions. These tests help ensure all functions
       related to JavaScript are working as intended.
    """

    def test_clear_file(self, selenium: webDriver):
        pass
