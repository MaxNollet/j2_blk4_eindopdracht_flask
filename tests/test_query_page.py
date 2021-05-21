from os import path
from typing import Union

from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.support import expected_conditions
from selenium.webdriver.support.select import Select
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


class TestDefaults:
    """A class which groups tests related to checking
       default values and settings for elements in the
       interface.
    """

    def test_defaults_query_generator(self, selenium: webDriver):
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

    def test_first_addition_no_spaces(self, selenium: webDriver):
        """Test the addition of a term to the query without spaces."""
        selenium.get("http://127.0.0.1:5000/")
        input_term = selenium.find_element_by_id("input_search_term")
        input_query = selenium.find_element_by_id("input_generated_query")
        button_add = selenium.find_element_by_id("button_add_item")
        term = "Pikanto"
        input_term.send_keys(term)
        selenium.execute_script("arguments[0].click();", button_add)
        assert input_term.get_attribute("value") == ""
        assert input_query.get_attribute("value") == f"{term}[ALL]"

    def test_first_addition_with_spaces(self, selenium: webDriver):
        """Test the addition of a term to the query with spaces."""
        selenium.get("http://127.0.0.1:5000/")
        input_term = selenium.find_element_by_id("input_search_term")
        input_query = selenium.find_element_by_id("input_generated_query")
        button_add = selenium.find_element_by_id("button_add_item")
        term = "Frikandel speciaal"
        input_term.send_keys(term)
        selenium.execute_script("arguments[0].click();", button_add)
        assert input_query.get_attribute("value") == f"{term}[ALL]"

    def test_second_addition_no_spaces(self, selenium: webDriver):
        """Test the addition of two terms with no spaces."""
        selenium.get("http://127.0.0.1:5000/")
        input_term = selenium.find_element_by_id("input_search_term")
        input_query = selenium.find_element_by_id("input_generated_query")
        button_add = selenium.find_element_by_id("button_add_item")
        term1 = "Kroket"
        term2 = "Frikandel"
        input_term.send_keys(term1)
        selenium.execute_script("arguments[0].click(0);", button_add)
        input_term.send_keys(term2)
        selenium.execute_script("arguments[0].click(0);", button_add)
        assert input_term.get_attribute("value") == ""
        assert input_query.get_attribute("value") == f"({term1}[ALL]) AND ({term2}[ALL])"

    def test_second_addition_with_spaces(self, selenium: webDriver):
        """Test the addition of two terms with spaces."""
        selenium.get("http://127.0.0.1:5000/")
        input_term = selenium.find_element_by_id("input_search_term")
        input_query = selenium.find_element_by_id("input_generated_query")
        button_add = selenium.find_element_by_id("button_add_item")
        term1 = "Frikandel speciaal"
        term2 = "Patatje oorlog"
        input_term.send_keys(term1)
        selenium.execute_script("arguments[0].click(0);", button_add)
        input_term.send_keys(term2)
        selenium.execute_script("arguments[0].click(0);", button_add)
        assert input_term.get_attribute("value") == ""
        assert input_query.get_attribute("value") == f"({term1}[ALL]) AND ({term2}[ALL])"

    def test_and_or_not(self, selenium: webDriver):
        """Test the addition of several items with different types of additions."""
        selenium.get("http://127.0.0.1:5000/")
        input_term = selenium.find_element_by_id("input_search_term")
        input_query = selenium.find_element_by_id("input_generated_query")
        input_add_type = Select(selenium.find_element_by_id("input_add_type"))
        term1 = "Viandel"
        term2 = "Bitterbal"
        term3 = "kipcorn"
        term4 = "Kaassoufflé"
        input_term.send_keys(term1, Keys.ENTER, term2, Keys.ENTER)
        input_add_type.select_by_value("OR")
        input_term.send_keys(term3, Keys.ENTER)
        input_add_type.select_by_value("NOT")
        input_term.send_keys(term4, Keys.ENTER)
        print(input_query.get_attribute("value"))
        assert input_term.get_attribute("value") == ""
        assert input_query.get_attribute(
            "value") == f"((({term1}[ALL]) AND ({term2}[ALL])) OR ({term3}[ALL])) NOT ({term4}[ALL])"

    def test_different_fields(self, selenium: webDriver):
        """Test the addition of multiple terms for different fields."""
        selenium.get("http://127.0.0.1:5000/")
        input_field = Select(selenium.find_element_by_id("input_field"))
        input_term = selenium.find_element_by_id("input_search_term")
        input_query = selenium.find_element_by_id("input_generated_query")
        input_add_type = Select(selenium.find_element_by_id("input_add_type"))
        fields = ("ALL", "BOOK", "ED", "PDAT", "JOUR", "SUBS")
        terms = ("Kipcorn", "Van-alles-wat", "Gehaktbal", "Gehaktbal speciaal", "Nasibal", "Bamischijf")
        additions = ("AND", "OR", "OR", "NOT", "NOT")
        input_field.select_by_value(fields[0])
        input_term.send_keys(terms[0], Keys.ENTER)
        for field, term, addition in zip(fields[1:], terms[1:], additions):
            input_field.select_by_value(field)
            input_add_type.select_by_value(addition)
            input_term.send_keys(term, Keys.ENTER)
        check = "(((((Kipcorn[ALL]) AND (Van-alles-wat[BOOK])) OR (Gehaktbal[ED])) OR (Gehaktbal speciaal[PDAT]))" \
                " NOT (Nasibal[JOUR])) NOT (Bamischijf[SUBS])"
        assert input_query.get_attribute("value") == check


class TestJavaScript:
    """A class which groups tests related to JavaScript-
       functions. These tests help ensure all functions
       related to JavaScript are working as intended.
    """

    def test_add_term_on_enter(self, selenium: webDriver):
        """Test the addition of a term to the query using the ENTER-key."""
        selenium.get("http://127.0.0.1:5000/")
        input_search_term = selenium.find_element_by_id("input_search_term")
        input_generated_query = selenium.find_element_by_id("input_generated_query")
        input_add_type = selenium.find_element_by_id("input_add_type")
        input_search_term.send_keys("Mexicano", Keys.ENTER)
        assert input_search_term.get_attribute("value") == ""
        assert input_generated_query.get_attribute("value") != ""
        assert input_add_type.is_enabled() is True

    def test_clear_file(self, selenium: webDriver):
        """Test the clear-button to remove the specified file."""
        selenium.get("http://127.0.0.1:5000/")
        input_load_symbols = selenium.find_element_by_id("input_load_symbols")
        button_clear_file = selenium.find_element_by_id("button_clear_file")
        file_path = path.abspath(path.join(path.dirname(__file__), "..", "requirements.txt"))
        input_load_symbols.send_keys(file_path)
        assert path.basename(input_load_symbols.get_attribute("value")) == path.basename(file_path)
        WebDriverWait(selenium, 10).until(
            expected_conditions.element_to_be_clickable((By.ID, "input_load_symbols"))
        )
        selenium.execute_script("arguments[0].click();", button_clear_file)
        assert input_load_symbols.get_attribute("value") == ""

    def test_open_new_tab(self, selenium: webDriver):
        """Test the open-in-new-tab-button."""
        selenium.get("http://127.0.0.1:5000/")
        check_new_tab = selenium.find_element_by_id("check_new_tab")
        input_form = selenium.find_element_by_id("input_form")
        assert check_new_tab.is_selected() is False
        assert input_form.get_attribute("target") == ""
        selenium.execute_script("arguments[0].click();", check_new_tab)
        assert check_new_tab.is_selected() is True
        assert input_form.get_attribute("target") == "_blank"
        selenium.execute_script("arguments[0].click();", check_new_tab)
        assert check_new_tab.is_selected() is False
        assert input_form.get_attribute("target") == ""


class TestEntireForm:
    """A class which groups methods required for testing the
       entire form at once.
    """

    def test_entire_form(self, selenium: webDriver):
        """Test the entire form."""
        selenium.maximize_window()
        selenium.get("http://127.0.0.1:5000/")
        selenium.execute_script("window.scrollTo(0,document.body.scrollHeight)")
        check_new_tab = selenium.find_element_by_id("check_new_tab")
        button_clear = selenium.find_element_by_id("button_clear")
        # Dry-run.
        self.build_query(selenium)
        self.specify_genes_exclude(selenium)
        self.specify_genes_include(selenium)
        self.optional_options(selenium)
        selenium.execute_script("arguments[0].click();", check_new_tab)
        selenium.execute_script("arguments[0].click();", button_clear)
        # Real run.
        query = self.build_query(selenium)
        self.assert_query_builder(selenium, query)
        self.specify_genes_exclude(selenium)
        self.assert_filename(selenium)
        symbols = self.specify_genes_include(selenium)
        self.assert_symbols(selenium, symbols)
        self.optional_options(selenium)

    @classmethod
    def build_query(cls, selenium: webDriver):
        """A method which build a query using the query
           builder in the form with various fields, terms
           and addition types.

        Input = selenium webdriver (WebDriver).
        """
        input_field = Select(selenium.find_element_by_id("input_field"))
        input_search_term = selenium.find_element_by_id("input_search_term")
        input_add_type = Select(selenium.find_element_by_id("input_add_type"))
        button_add_item = selenium.find_element_by_id("button_add_item")
        fields = ("AFFL", "ALL", "AUTH", "DSO", "CRDT", "EID", "FILT", "ISBN", "MESH", "PTYP", "WORD", "UID")
        terms = ("Inkscape", "GIMP", "Autodesk 3ds Max", "Autodesk Maya", "Autodesk SketchBook", "HitFilm",
                 "SceneBuilder", "Bootstrap Studio", "Jetbrains IntelliJ IDEA", "Jetbrains PyCharm",
                 "Microsoft Visual Studio", "Windows Terminal")
        additions = ("AND", "AND", "OR", "AND", "NOT", "NOT", "AND", "OR", "OR", "AND", "NOT")
        input_field.select_by_value(fields[0])
        input_search_term.send_keys(terms[0])
        selenium.execute_script("arguments[0].click();", button_add_item)
        for field, term, addition in zip(fields[1:], terms[1:], additions):
            input_field.select_by_value(field)
            input_add_type.select_by_value(addition)
            input_search_term.send_keys(term)
            selenium.execute_script("arguments[0].click();", button_add_item)
        query = "(((((((((((Inkscape[AFFL]) AND (GIMP[ALL])) AND (Autodesk 3ds Max[AUTH])) OR (Autodesk Maya[DSO]))" \
                " AND (Autodesk SketchBook[CRDT])) NOT (HitFilm[EID])) NOT (SceneBuilder[FILT])) AND" \
                " (Bootstrap Studio[ISBN])) OR (Jetbrains IntelliJ IDEA[MESH])) OR (Jetbrains PyCharm[PTYP])) AND" \
                " (Microsoft Visual Studio[WORD])) NOT (Windows Terminal[UID])"
        return query

    @classmethod
    def specify_genes_include(cls, selenium: webDriver):
        """A method which enters various genes to be included
           into the search.

        Input = selenium webdriver (WebDriver).
        """
        radio_include_symbols = selenium.find_element_by_id("radio_include_symbols")
        radio_exclude_symbols = selenium.find_element_by_id("radio_exclude_symbols")
        input_symbols = selenium.find_element_by_id("input_symbols")
        symbols = ("ABCA12", "KPPP1", "A4GALT", "EST140535", "EST349056", "ABCC1", "MRP8", "ABCG5", "CALJA", "CEPU-1",
                   "POMP", "RNF168", "MYMY2", "ORP1", "DBA19", "RPS15A")
        selenium.execute_script("arguments[0].click();", radio_include_symbols)
        selenium.execute_script("arguments[0].click();", radio_exclude_symbols)
        selenium.execute_script("arguments[0].click();", radio_include_symbols)
        for symbol in symbols[:-1]:
            input_symbols.send_keys(symbol, ", ")
        input_symbols.send_keys(symbols[-1])
        return ", ".join(symbols)

    @classmethod
    def specify_genes_exclude(cls, selenium: webDriver):
        """A method which uploads a file to exclude
           from the search.

        Input = selenium webdriver (WebDriver).
        """
        radio_include_symbols = selenium.find_element_by_id("radio_include_symbols")
        radio_exclude_symbols = selenium.find_element_by_id("radio_exclude_symbols")
        input_load_symbols = selenium.find_element_by_id("input_load_symbols")
        button_clear_file = selenium.find_element_by_id("button_clear_file")
        selenium.execute_script("arguments[0].click();", radio_exclude_symbols)
        selenium.execute_script("arguments[0].click();", radio_include_symbols)
        selenium.execute_script("arguments[0].click();", radio_exclude_symbols)
        file = path.abspath(__file__)
        input_load_symbols.send_keys(file)
        WebDriverWait(selenium, 10).until(
            expected_conditions.text_to_be_present_in_element_value((By.ID, "input_load_symbols"),
                                                                    path.basename(file))
        )
        selenium.execute_script("arguments[0].click()", button_clear_file)

    @classmethod
    def optional_options(cls, selenium: webDriver):
        """A method which sets optional options.
        Input = selenium webdriver (WebDriver).
        """
        input_date_after = selenium.find_element_by_id("input_date_after")
        input_date_before = selenium.find_element_by_id("input_date_before")
        date1 = "24052000"
        date2 = "20052021"
        input_date_after.send_keys(date1)
        input_date_before.send_keys(date2)
        return [date1, date2]

    @classmethod
    def assert_query_builder(cls, selenium: webDriver, query):
        """Assert that every element in the query builder is ok."""
        assert selenium.find_element_by_id("input_search_term").is_enabled() is True
        assert selenium.find_element_by_id("input_add_type").is_enabled() is True
        assert selenium.find_element_by_id("input_generated_query").get_attribute("value") == query

    @classmethod
    def assert_symbols(cls, selenium: webDriver, symbols):
        """Assert that symbols are inserted and that the right
           options are selected.
        """
        radio_include_symbols = selenium.find_element_by_id("radio_include_symbols")
        radio_exclude_symbols = selenium.find_element_by_id("radio_exclude_symbols")
        input_symbols = selenium.find_element_by_id("input_symbols")
        input_load_symbols = selenium.find_element_by_id("input_load_symbols")
        assert radio_include_symbols.is_selected() is True
        assert radio_exclude_symbols.is_selected() is False
        assert input_symbols.get_attribute("value") == symbols
        assert input_load_symbols.get_attribute("value") == ""

    @classmethod
    def assert_filename(cls, selenium: webDriver):
        """Assert that the file-upload is ok and that the right
           option are selected.
        """
        radio_include_symbols = selenium.find_element_by_id("radio_include_symbols")
        radio_exclude_symbols = selenium.find_element_by_id("radio_exclude_symbols")
        input_symbols = selenium.find_element_by_id("input_symbols")
        input_load_symbols = selenium.find_element_by_id("input_load_symbols")
        assert radio_include_symbols.is_selected() is False
        assert radio_exclude_symbols.is_selected() is True
        assert input_symbols.get_attribute("value") == ""
        assert input_load_symbols.get_attribute("value") == ""
