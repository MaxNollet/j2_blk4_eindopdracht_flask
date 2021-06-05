/**
 * A function which adds several event triggers to
 * elements in the interface.
 *
 * @author Max Nollet
 * */
window.onload = function () {
    // Listener for input search term.
    document.getElementById("input_search_term")
        .addEventListener("keypress", function (event) {
            if (event.key === "Enter") {
                event.preventDefault();
                document.getElementById("button_add_item").click();
            }
        });
    // Listener for button to add the search term to the query.
    document.getElementById("button_add_item")
        .addEventListener("click", QueryBuilder);
    // Listener for changes in the query-field.
    document.getElementById("input_generated_query")
        .addEventListener("input", OnQueryChange);
    OnQueryChange();
    // Listener for button to clear the chosen file.
    document.getElementById("button_clear_file")
        .addEventListener("click", function () {
            ClearFile("input_load_symbols")
        });
    // Listener for checkbox to open results in a new tab.
    document.getElementById("check_new_tab")
        .addEventListener("click", function () {
            OpenNewTab(this.id, "input_form");
        });
    // Listener for reset of the form to reset open in new tab.
    document.getElementById("input_form")
        .addEventListener("reset", function () {
            OpenNewTab(this.id, "input_form");
        });
    document.getElementById("input_form")
        .addEventListener("submit", SubmitQuery);
}

/**
 * A function which retrieves values from the interface
 * and uses these values to update the interface with a
 * new query.
 *
 * @author Max Nollet
 * */
function QueryBuilder() {
    const element_field = document.getElementById("input_field");
    const element_term = document.getElementById("input_search_term");
    const element_type = document.getElementById("input_add_type");
    const element_query = document.getElementById("input_generated_query");
    const field = element_field.options[element_field.selectedIndex].value;
    const term = element_term.value.trim();
    const type = element_type.options[element_type.selectedIndex].value;
    const query = element_query.value.trim();
    if (term !== "") {
        element_query.value = Concatenate(field, term, type, query);
        element_term.value = "";
        OnQueryChange();
    }
}

/**
 * A function which concatenates values to build a new
 * query.
 *
 * @author Max Nollet
 * @param field Field to be searched in scientific articles.
 * @param term Search term entered by the user to search scientific articles.
 * @param type Type of addition to the query (AND, OR, NOT).
 * @param query Base query where new values are appended.
 * @return query New concatenated query.
 * */
function Concatenate(field, term, type, query) {
    if (term !== "") {
        if (query === "") {
            query = `${term}[${field}]`;
        } else {
            query = `(${query}) ${type} (${term}[${field}])`;
        }
    }
    return query;
}

/**
 * A function which enables/disables the dropdown
 * button containing different types for adding a
 * value to the query.
 *
 * @author Max Nollet
 * */
function OnQueryChange() {
    const element_query = document.getElementById("input_generated_query");
    const element_add = document.getElementById("input_add_type");
    element_add.disabled = element_query.value.trim() === "";
}

/***/
function SubmitQuery(event) {
    event.preventDefault();
    RemoveChilds("alert_box");
    const original_text = document.getElementById("button_submit").innerText;
    ToggleDisableSubmitButton(true, "button_submit", original_text);

    const form_element = document.getElementById("input_form");

    const target_url = form_element.getAttribute("action");
    const request = new XMLHttpRequest();
    request.open("POST", target_url, true);
    // request.setRequestHeader("Content-Type", "application/x-www-form-urlencoded");
    request.onreadystatechange = function () {
        ResponseHandler(this, "alert_box");
        ToggleDisableSubmitButton(false, "button_submit", "Search genes");
    };
    request.send(new FormData(form_element));
}
