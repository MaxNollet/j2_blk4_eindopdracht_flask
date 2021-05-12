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
        const new_query = Concatenate(field, term, type, query);
        if (new_query !== "") {
            element_query.value = new_query;
            element_term.value = "";
        }
    }
}

function Concatenate(field, term, type, query) {
    if (term !== "") {
        if (query === "") {
            return `${term}[${field}]`;
        } else {
            return `(${query}) ${type} (${term}[${field}])`;
        }
    } else return "";
}
