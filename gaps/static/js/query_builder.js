/**
 * A function which adds several event triggers to
 * elements in the interface.
 *
 * @author Max Nollet
 * */
window.onload = function () {
    document.getElementById("input_search_term")
        .addEventListener("keypress", function (event) {
            if (event.key === "Enter") {
                event.preventDefault();
                document.getElementById("button_add_item").click();
            }
        });
    document.getElementById("button_add_item")
        .addEventListener("click", QueryBuilder);
    document.getElementById("input_generated_query")
        .addEventListener("change", OnQueryChange);
    OnQueryChange();
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

function OnQueryChange() {
    const element_query = document.getElementById("input_generated_query");
    const element_add = document.getElementById("input_add_type");
    element_add.disabled = element_query.value.trim() === "";
}
