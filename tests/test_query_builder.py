import ntpath
from os import path
from typing import Union

from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.remote.webelement import WebElement
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


class ElementSelector:
    """A class that groups methods which select various
       elements on the query_builder-page.

    base_url = URL of the page containing the elements (str).
    """
    base_url = "http://127.0.0.1:5000/query_builder"

    @staticmethod
    def select_input_field(selenium: webDriver):
        """Select element 'input_field'."""
        return selenium.find_element_by_id("input_field")

    @staticmethod
    def select_input_search_term(selenium: webDriver):
        """Select element 'input_search_term'."""
        return selenium.find_element_by_id("input_search_term")

    @staticmethod
    def select_input_add_type(selenium: webDriver):
        """Select element 'input_add_type'."""
        return selenium.find_element_by_id("input_add_type")

    @staticmethod
    def select_button_add_item(selenium: webDriver):
        """Select element 'button_add_item'."""
        return selenium.find_element_by_id("button_add_item")

    @staticmethod
    def select_input_generated_query(selenium: webDriver):
        """Select element 'input_generated_query'."""
        return selenium.find_element_by_id("input_generated_query")

    @staticmethod
    def select_radio_include_symbols(selenium: webDriver):
        """Select element 'radio_include_symbols'."""
        return selenium.find_element_by_id("radio_include_symbols")

    @staticmethod
    def select_radio_exclude_symbols(selenium: webDriver):
        """Select element 'radio_exclude_symbols'."""
        return selenium.find_element_by_id("radio_exclude_symbols")

    @staticmethod
    def select_input_symbols(selenium: webDriver):
        """Select element 'input_symbols'."""
        return selenium.find_element_by_id("input_symbols")

    @staticmethod
    def select_input_load_symbols(selenium: webDriver):
        """Select element 'input_load_symbols'."""
        return selenium.find_element_by_id("input_load_symbols")

    @staticmethod
    def select_button_clear_file(selenium: webDriver):
        """Select element 'button_clear_file'."""
        return selenium.find_element_by_id("button_clear_file")

    @staticmethod
    def select_input_date_after(selenium: webDriver):
        """Select element 'input_date_after'."""
        return selenium.find_element_by_id("input_date_after")

    @staticmethod
    def select_input_date_before(selenium: webDriver):
        """Select element 'input_date_before'."""
        return selenium.find_element_by_id("input_date_before")

    @staticmethod
    def select_check_new_tab(selenium: webDriver):
        """Select element 'check_new_tab'."""
        return selenium.find_element_by_id("check_new_tab")

    @staticmethod
    def select_button_submit(selenium: webDriver):
        """Select element 'button_submit'."""
        return selenium.find_element_by_id("button_submit")

    @staticmethod
    def select_button_clear(selenium: webDriver):
        """Select element 'button_clear'."""
        return selenium.find_element_by_id("button_clear")

    @staticmethod
    def select_input_form(selenium: webDriver):
        """Select element 'input_form'."""
        return selenium.find_element_by_id("input_form")


class HelperFunctions:
    """A class that groups functions which help execute
       something to prevent typing the same stuff over
       and over again.
    """

    @staticmethod
    def click_element(selenium: webDriver, element: WebElement) -> None:
        """A method which clicks an element using a JavaScript-
           function.

        :param selenium Selenium webdriver (WebDriver).
        :param element Element to be clicked on (WebElement).
        """
        selenium.execute_script("arguments[0].click();", element)

    @classmethod
    def build_query(cls, selenium: webDriver, fields: tuple, terms: tuple, additions: tuple = None) -> None:
        """A method which helps to build a query from lists of
           values.

        :param selenium Selenium webdriver (WebDriver).
        :param fields Fields to search for (tuple).
        :param terms Terms to search (tuple).
        :param additions Types of additions to use (tuple).
        """
        input_field = Select(ElementSelector.select_input_field(selenium))
        input_search_term = ElementSelector.select_input_search_term(selenium)
        input_add_type = Select(ElementSelector.select_input_add_type(selenium))
        input_field.select_by_value(fields[0])
        input_search_term.send_keys(terms[0], Keys.ENTER)
        if len(fields) > 1 and len(terms) > 1 and additions is not None:
            for field, term, addition in zip(fields[1:], terms[1:], additions):
                input_field.select_by_value(field)
                input_add_type.select_by_value(addition)
                input_search_term.send_keys(term, Keys.ENTER)


class TestDefaultsQueryGenerator(ElementSelector):
    """A class which groups tests that check the defaults
       for the query builder in step 1 in the interface.
    """

    # input_field
    def test_input_field_enabled(self, selenium: webDriver):
        """Test input_field if enabled."""
        selenium.get(self.base_url)
        assert self.select_input_field(selenium).is_enabled() is True

    def test_input_field_texts(self, selenium: webDriver):
        """Test input_field if all text-values are correct."""
        selenium.get(self.base_url)
        texts = ("All Fields", "UID", "Filter", "Title", "Text Word", "MeSH Terms", "MeSH Major Topic", "Author",
                 "Journal", "Affiliation", "EC/RN Number", "Supplementary Concept", "Date - Publication",
                 "Date - Entrez", "Volume", "Pagination", "Publication Type", "Language", "Issue", "MeSH Subheading",
                 "Secondary Source ID", "Date - MeSH", "Title/Abstract", "Other Term", "Investigator",
                 "Author - Corporate", "Place of Publication", "Pharmacological Action", "Grant Number",
                 "Date - Modification", "Date - Completion", "Publisher ID", "Author - First", "Author - Full",
                 "Investigator - Full", "Transliterated Title", "Author - Last", "Print Publication Date",
                 "Electronic Publication Date", "Location ID", "Date - Create", "Book", "Editor", "ISBN", "Publisher",
                 "Author Cluster ID", "Extended PMID", "DSO", "Author - Identifier", "Subject - Personal Name",
                 "Conflict of Interest Statements")
        assert all(option.text in texts for option in Select(self.select_input_field(selenium)).options)

    def test_input_field_values(self, selenium: webDriver):
        """Test input_field if all values are correct"""
        selenium.get(self.base_url)
        values = ("ALL", "UID", "FILT", "TITL", "WORD", "MESH", "MAJR", "AUTH", "JOUR", "AFFL", "ECNO", "SUBS", "PDAT",
                  "EDAT", "VOL", "PAGE", "PTYP", "LANG", "ISS", "SUBH", "SI", "MHDA", "TIAB", "OTRM", "INVR", "COLN",
                  "CNTY", "PAPX", "GRNT", "MDAT", "CDAT", "PID", "FAUT", "FULL", "FINV", "TT", "LAUT", "PPDT", "EPDT",
                  "LID", "CRDT", "BOOK", "ED", "ISBN", "PUBN", "AUCL", "EID", "DSO", "AUID", "PS", "COIS")
        assert all(option.get_attribute("value") in values for option
                   in Select(self.select_input_field(selenium)).options)

    def test_input_field_default(self, selenium: webDriver):
        """Test input_field default value."""
        selenium.get(self.base_url)
        assert Select(self.select_input_field(selenium)).first_selected_option.text == "All Fields"

    # input_search_term
    def test_input_search_term_enabled(self, selenium: webDriver):
        """Test input_search_term if enabled."""
        selenium.get(self.base_url)
        assert self.select_input_search_term(selenium).is_enabled() is True

    def test_input_search_term_default(self, selenium: webDriver):
        """Test input_search_term default value."""
        selenium.get(self.base_url)
        assert self.select_input_search_term(selenium).get_attribute("value") == ""

    # input_add_type
    def test_input_add_type_enabled(self, selenium: webDriver):
        """Test input_add_type if enabled."""
        selenium.get(self.base_url)
        assert self.select_input_add_type(selenium).is_enabled() is False

    def test_input_add_type_texts(self, selenium: webDriver):
        """Test input_add_type if all text-values are correct."""
        selenium.get(self.base_url)
        texts = ("AND", "OR", "NOT")
        assert all(option.text in texts for option in Select(self.select_input_add_type(selenium)).options)

    def test_input_add_type_values(self, selenium: webDriver):
        """Test input_add_type if all values are correct."""
        selenium.get(self.base_url)
        options = ("AND", "OR", "NOT")
        assert all(option.get_attribute("value") in options for option
                   in Select(self.select_input_add_type(selenium)).options)

    def test_input_add_type_default(self, selenium: webDriver):
        """Test input_add_type default value."""
        selenium.get(self.base_url)
        assert Select(self.select_input_add_type(selenium)).first_selected_option.text == "AND"

    # button_add_item
    def test_button_add_item_enabled(self, selenium: webDriver):
        """Test button_add_item if enabled."""
        selenium.get(self.base_url)
        assert self.select_button_add_item(selenium).is_enabled() is True

    def test_button_add_item_text(self, selenium: webDriver):
        """Test button_add_item if text is correct"""
        selenium.get(self.base_url)
        assert self.select_button_add_item(selenium).text == "Add"

    # input_generated_query
    def test_input_generated_query_enabled(self, selenium: webDriver):
        """Test input_generated_query if enabled."""
        selenium.get(self.base_url)
        assert self.select_input_generated_query(selenium).is_enabled() is True

    def test_input_generated_query_default(self, selenium: webDriver):
        """Test input_generated_query default value."""
        selenium.get(self.base_url)
        assert self.select_input_generated_query(selenium).get_attribute("value") == ""


class TestDefaultsGeneSymbols(ElementSelector):
    """A class which groups tests that check the defaults
       for the gene symbols section in step 2 in the interface.
    """

    # radio_include_symbols
    def test_radio_include_symbols_enabled(self, selenium: webDriver):
        """Test radio_include_symbols if enabled."""
        selenium.get(self.base_url)
        assert self.select_radio_include_symbols(selenium).is_enabled() is True

    def test_radio_include_symbols_selected(self, selenium: webDriver):
        """Test radio_include_symbols if selected."""
        selenium.get(self.base_url)
        assert self.select_radio_include_symbols(selenium).is_selected() is True

    def test_radio_include_symbols_value(self, selenium: webDriver):
        """Test radio_include_symbols value."""
        selenium.get(self.base_url)
        assert self.select_radio_include_symbols(selenium).get_attribute("value") == "true"

    # radio_exclude_symbols
    def test_radio_exclude_symbols_enabled(self, selenium: webDriver):
        """Test radio_exclude_symbols if enabled."""
        selenium.get(self.base_url)
        assert self.select_radio_exclude_symbols(selenium).is_enabled() is True

    def test_radio_exclude_symbols_selected(self, selenium: webDriver):
        """Test radio_exclude_symbols if selected."""
        selenium.get(self.base_url)
        assert self.select_radio_exclude_symbols(selenium).is_selected() is False

    def test_radio_exclude_symbols_value(self, selenium: webDriver):
        """Test radio_exclude_symbols value."""
        selenium.get(self.base_url)
        assert self.select_radio_exclude_symbols(selenium).get_attribute("value") == "false"

    # input_symbols
    def test_input_symbols_enabled(self, selenium: webDriver):
        """Test input_symbols is enabled."""
        selenium.get(self.base_url)
        assert self.select_input_symbols(selenium).is_enabled() is True

    def test_input_symbols_default(self, selenium: webDriver):
        """Test input_symbols default value."""
        selenium.get(self.base_url)
        assert self.select_input_symbols(selenium).get_attribute("value") == ""

    # input_load_symbols
    def test_input_load_symbols_enabled(self, selenium: webDriver):
        """Test input_load_symbols if enabled."""
        selenium.get(self.base_url)
        assert self.select_input_load_symbols(selenium).is_enabled() is True

    def test_input_load_symbols_default(self, selenium: webDriver):
        """Test input_load_symbols default value."""
        selenium.get(self.base_url)
        assert self.select_input_load_symbols(selenium).get_attribute("value") == ""

    # button_clear_file
    def test_button_clear_file_enabled(self, selenium: webDriver):
        """Test button_clear_file is enabled."""
        selenium.get(self.base_url)
        assert self.select_button_clear_file(selenium).is_enabled() is True

    def test_button_clear_file_text(self, selenium: webDriver):
        """Test button_clear_file if text is correct."""
        selenium.get(self.base_url)
        assert self.select_button_clear_file(selenium).text == "Clear"


class TestDefaultsOptionalOptions(ElementSelector):
    """A class which groups tests that check the defaults
       for the optional options in step 3 in the interface.
    """

    # input_date_after
    def test_input_date_after_enabled(self, selenium: webDriver):
        """Test input_date_after if enabled."""
        selenium.get(self.base_url)
        assert self.select_input_date_after(selenium).is_enabled() is True

    def test_input_date_after_default(self, selenium: webDriver):
        """Test input_date_after default value."""
        selenium.get(self.base_url)
        assert self.select_input_date_after(selenium).get_attribute("value") == ""

    # input_date_before
    def test_input_date_before_enabled(self, selenium: webDriver):
        """Test input_date_before if enabled."""
        selenium.get(self.base_url)
        assert self.select_input_date_before(selenium).is_enabled() is True

    def test_input_date_before_default(self, selenium: webDriver):
        """Test input_date_before default value."""
        selenium.get(self.base_url)
        assert self.select_input_date_before(selenium).get_attribute("value") == ""

    # check_new_tab
    def test_check_new_tab_enabled(self, selenium: webDriver):
        """Test check_new_tab if enabled."""
        selenium.get(self.base_url)
        assert self.select_check_new_tab(selenium).is_enabled() is True

    def test_check_new_tab_selected(self, selenium: webDriver):
        """Test check_new_tab if selected."""
        selenium.get(self.base_url)
        assert self.select_check_new_tab(selenium).is_selected() is False

    # button_submit
    def test_button_submit_enabled(self, selenium: webDriver):
        """Test button_submit if enabled."""
        selenium.get(self.base_url)
        assert self.select_button_submit(selenium).is_enabled() is True

    def test_button_submit_type(self, selenium: webDriver):
        """Test button_submit if correct type."""
        selenium.get(self.base_url)
        assert self.select_button_submit(selenium).get_attribute("type") == "submit"

    # button_clear
    def test_button_clear_enabled(self, selenium: webDriver):
        """Test button_clear is enabled."""
        selenium.get(self.base_url)
        assert self.select_button_clear(selenium).is_enabled() is True

    def test_button_clear_type(self, selenium: webDriver):
        """Test button_clear if correct type."""
        selenium.get(self.base_url)
        assert self.select_button_clear(selenium).get_attribute("type") == "reset"


class TestQueryBuilder(ElementSelector, HelperFunctions):
    """A class which groups tests related to the query
       builder. These tests are in place to ensure the
       query builder keeps working as intended.
    """

    # first addition.
    def test_one_addition_on_enter(self, selenium: webDriver):
        """Test one addition to a query by using the 'ENTER'-key."""
        selenium.get(self.base_url)
        input_search_term = self.select_input_search_term(selenium)
        term = "Pikanto"
        input_search_term.send_keys(term, Keys.ENTER)
        assert input_search_term.get_attribute("value") == ""

    def test_one_addition_on_click(self, selenium: webDriver):
        """Test one addition to a query by clicking on the 'Add'-button."""
        selenium.get(self.base_url)
        input_search_term = self.select_input_search_term(selenium)
        button_add_item = self.select_button_add_item(selenium)
        term = "Pikanto"
        input_search_term.send_keys(term)
        self.click_element(selenium, button_add_item)
        assert input_search_term.get_attribute("value") == ""

    def test_one_addition_result(self, selenium: webDriver):
        """Test one addition to a query and check the result."""
        selenium.get(self.base_url)
        input_search_term = self.select_input_search_term(selenium)
        input_generated_query = self.select_input_generated_query(selenium)
        term = "Pikanto"
        input_search_term.send_keys(term, Keys.ENTER)
        assert input_generated_query.get_attribute("value") == f"{term}[ALL]"

    def test_one_addition_enable_addition(self, selenium: webDriver):
        """Test one addition to a query and check if the addition types are enabled."""
        selenium.get(self.base_url)
        input_search_term = self.select_input_search_term(selenium)
        input_add_type = self.select_input_add_type(selenium)
        term = "Pikanto"
        input_search_term.send_keys(term, Keys.ENTER)
        assert input_add_type.is_enabled() is True

    # Second addition.
    def test_two_addition_and(self, selenium: webDriver):
        """Test two additions of two terms to the query with an 'And'-addition."""
        selenium.get(self.base_url)
        fields = ("ALL", "ALL")
        terms = ("Kroket", "Frikandel")
        additions = ("AND",)
        self.build_query(selenium, fields, terms, additions)
        input_generated_query = self.select_input_generated_query(selenium)
        assert input_generated_query.get_attribute("value") == f"({terms[0]}[ALL]) AND ({terms[1]}[ALL])"

    def test_second_addition_with_spaces(self, selenium: webDriver):
        """Test two additions of two elements to the query with and 'Or'-addition."""
        selenium.get(self.base_url)
        fields = ("ALL", "ALL")
        terms = ("Frikandel speciaal", "Patatje oorlog")
        additions = ("OR",)
        self.build_query(selenium, fields, terms, additions)
        input_generated_query = self.select_input_generated_query(selenium)
        assert input_generated_query.get_attribute("value") == f"({terms[0]}[ALL]) OR ({terms[1]}[ALL])"

    # Miscellaneous tests.
    def test_and_or_not(self, selenium: webDriver):
        """Test the addition of several items with different types of additions."""
        selenium.get(self.base_url)
        fields = ("ALL", "ALL", "ALL", "ALL")
        terms = ("Viandel", "Bitterbal", "kipcorn", "Kaassoufflé")
        additions = ("AND", "OR", "NOT")
        self.build_query(selenium, fields, terms, additions)
        input_generated_query = self.select_input_generated_query(selenium)
        assert input_generated_query.get_attribute("value") == f"((({terms[0]}[ALL]) AND ({terms[1]}[ALL])) OR " \
                                                               f"({terms[2]}[ALL])) NOT ({terms[3]}[ALL])"

    def test_different_fields_terms_additions(self, selenium: webDriver):
        """Test the addition of multiple terms for different fields."""
        selenium.get(self.base_url)
        fields = ("ALL", "BOOK", "ED", "PDAT", "JOUR", "SUBS")
        terms = ("Kipcorn", "Van-alles-wat", "Gehaktbal", "Gehaktbal speciaal", "Nasibal", "Bamischijf")
        additions = ("AND", "OR", "OR", "NOT", "NOT")
        check = "(((((Kipcorn[ALL]) AND (Van-alles-wat[BOOK])) OR (Gehaktbal[ED])) OR (Gehaktbal speciaal[PDAT]))" \
                " NOT (Nasibal[JOUR])) NOT (Bamischijf[SUBS])"
        self.build_query(selenium, fields, terms, additions)
        input_generated_query = self.select_input_generated_query(selenium)
        assert input_generated_query.get_attribute("value") == check

    def test_big_query(self, selenium: webDriver):
        """Test the addition of multiple fields, terms and types of additions."""
        selenium.get(self.base_url)
        fields = ("AFFL", "ALL", "AUTH", "DSO", "CRDT", "EID", "FILT", "ISBN", "MESH", "PTYP", "WORD", "UID")
        terms = ("Inkscape", "GIMP", "Autodesk 3ds Max", "Autodesk Maya", "Autodesk SketchBook", "HitFilm",
                 "SceneBuilder", "Bootstrap Studio", "Jetbrains IntelliJ IDEA", "Jetbrains PyCharm",
                 "Microsoft Visual Studio", "Windows Terminal")
        additions = ("AND", "AND", "OR", "AND", "NOT", "NOT", "AND", "OR", "OR", "AND", "NOT")
        check = "(((((((((((Inkscape[AFFL]) AND (GIMP[ALL])) AND (Autodesk 3ds Max[AUTH])) OR (Autodesk Maya[DSO]))" \
                " AND (Autodesk SketchBook[CRDT])) NOT (HitFilm[EID])) NOT (SceneBuilder[FILT])) AND" \
                " (Bootstrap Studio[ISBN])) OR (Jetbrains IntelliJ IDEA[MESH])) OR (Jetbrains PyCharm[PTYP])) AND" \
                " (Microsoft Visual Studio[WORD])) NOT (Windows Terminal[UID])"
        self.build_query(selenium, fields, terms, additions)
        input_generated_query = self.select_input_generated_query(selenium)
        assert input_generated_query.get_attribute("value") == check


class TestJavaScript(ElementSelector, HelperFunctions):
    """A class which groups tests related to JavaScript-
       functions. These tests help ensure all functions
       related to JavaScript are working as intended.
    """

    def test_file_upload(self, selenium: webDriver):
        """Test for uploading files."""
        selenium.get(self.base_url)
        input_load_symbols = selenium.find_element_by_id("input_load_symbols")
        file_path = path.abspath(path.join(path.dirname(__file__), "..", "requirements.txt"))
        input_load_symbols.send_keys(file_path)
        assert ntpath.basename(input_load_symbols.get_attribute("value")) == ntpath.basename(file_path)

    def test_file_clear(self, selenium: webDriver):
        """Test the clear-button to remove the specified file."""
        selenium.get(self.base_url)
        input_load_symbols = selenium.find_element_by_id("input_load_symbols")
        button_clear_file = selenium.find_element_by_id("button_clear_file")
        file_path = path.abspath(path.join(path.dirname(__file__), "..", "requirements.txt"))
        input_load_symbols.send_keys(file_path)
        self.click_element(selenium, button_clear_file)
        assert input_load_symbols.get_attribute("value") == ""

    def test_open_new_tab_select(self, selenium: webDriver):
        """Test the open-in-new-tab-button when selected."""
        selenium.get(self.base_url)
        check_new_tab = self.select_check_new_tab(selenium)
        input_form = self.select_input_form(selenium)
        self.click_element(selenium, check_new_tab)
        assert input_form.get_attribute("target") == "_blank"

    def test_open_new_tab_deselect(self, selenium: webDriver):
        """Test the open-in-new-tab-button when deselected."""
        selenium.get(self.base_url)
        check_new_tab = self.select_check_new_tab(selenium)
        input_form = self.select_input_form(selenium)
        self.click_element(selenium, check_new_tab)
        self.click_element(selenium, check_new_tab)
        assert input_form.get_attribute("target") == ""

    def test_open_new_tab_reset(self, selenium: webDriver):
        """Test the open-in-new-tab-button with a complete form reset."""
        selenium.get(self.base_url)
        check_new_tab = self.select_check_new_tab(selenium)
        input_form = self.select_input_form(selenium)
        button_clear = self.select_button_clear(selenium)
        self.click_element(selenium, check_new_tab)
        self.click_element(selenium, button_clear)
        assert input_form.get_attribute("target") == ""


class TestEntireForm(ElementSelector, HelperFunctions):
    """A class which groups methods for testing several
       elements in the form.
    """

    def test_radio_group_include_selected(self, selenium: webDriver):
        """Test the radio-button 'include' when selected."""
        selenium.get(self.base_url)
        radio_include_symbols = self.select_radio_include_symbols(selenium)
        radio_exclude_symbols = self.select_radio_exclude_symbols(selenium)
        self.click_element(selenium, radio_include_symbols)
        self.click_element(selenium, radio_exclude_symbols)
        self.click_element(selenium, radio_include_symbols)
        assert radio_include_symbols.is_selected() is True

    def test_radio_group_include_deselect(self, selenium: webDriver):
        """Test the radio-button 'include' when not selected."""
        selenium.get(self.base_url)
        radio_include_symbols = self.select_radio_include_symbols(selenium)
        radio_exclude_symbols = self.select_radio_exclude_symbols(selenium)
        self.click_element(selenium, radio_include_symbols)
        self.click_element(selenium, radio_exclude_symbols)
        self.click_element(selenium, radio_include_symbols)
        assert radio_exclude_symbols.is_selected() is False

    def test_radio_group_exclude_selected(self, selenium: webDriver):
        """Test the radio-button 'exclude' when selected."""
        selenium.get(self.base_url)
        radio_include_symbols = self.select_radio_include_symbols(selenium)
        radio_exclude_symbols = self.select_radio_exclude_symbols(selenium)
        self.click_element(selenium, radio_exclude_symbols)
        self.click_element(selenium, radio_include_symbols)
        self.click_element(selenium, radio_exclude_symbols)
        assert radio_exclude_symbols.is_selected() is True

    def test_radio_group_exclude_deselected(self, selenium: webDriver):
        """Test the radio-button 'exclude' when not selected."""
        selenium.get(self.base_url)
        radio_include_symbols = self.select_radio_include_symbols(selenium)
        radio_exclude_symbols = self.select_radio_exclude_symbols(selenium)
        self.click_element(selenium, radio_exclude_symbols)
        self.click_element(selenium, radio_include_symbols)
        self.click_element(selenium, radio_exclude_symbols)
        assert radio_include_symbols.is_selected() is False

    def test_include_genes(self, selenium: webDriver):
        """Test the input of gene symbols."""
        selenium.get(self.base_url)
        radio_include_symbols = self.select_radio_include_symbols(selenium)
        radio_exclude_symbols = self.select_radio_exclude_symbols(selenium)
        input_symbols = self.select_input_symbols(selenium)
        symbols = ("ABCA12", "KPPP1", "A4GALT", "EST140535", "EST349056", "ABCC1", "MRP8", "ABCG5", "CALJA", "CEPU-1",
                   "POMP", "RNF168", "MYMY2", "ORP1", "DBA19", "RPS15A")
        self.click_element(selenium, radio_include_symbols)
        self.click_element(selenium, radio_exclude_symbols)
        self.click_element(selenium, radio_include_symbols)
        counter = 0
        for symbol in symbols:
            counter += 1
            input_symbols.send_keys(symbol)
            if counter < len(symbols):
                input_symbols.send_keys(", ")
        assert input_symbols.get_attribute("value") == ", ".join(symbols)
