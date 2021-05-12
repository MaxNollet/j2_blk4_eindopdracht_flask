/**
 * A function which adds several event triggers to
 * elements in the interface.
 *
 * @author Max Nollet
 * */
window.onload = function () {
    document.getElementById("query_search_term")
        .addEventListener("keypress", function (event) {
            if (event.key === "Enter") {
                event.preventDefault();
                document.getElementById("add_search_item").click();
            }
        });
    document.getElementById("add_search_item")
        .addEventListener("click", QueryBuilder);
}

/**
 * A function which retrieves values from the interface
 * and uses these values to update the interface with a
 * new query.
 *
 * @author Max Nollet
 * */
function QueryBuilder() {
    const element_field = document.getElementById("query_field");
    const element_term = document.getElementById("query_search_term");
    const element_type = document.getElementById("query_add_type");
    const element_query = document.getElementById("query_generated");
    const field = element_field.options[element_field.selectedIndex].value;
    const term = element_term.value.trim();
    const type = element_type.options[element_type.selectedIndex].value;
    const query = element_query.value.trim();
    if (term !== "") {
        element_query.value = Concatenate(field, term, type, query);
        element_term.value = "";
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
