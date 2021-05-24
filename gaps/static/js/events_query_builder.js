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
}
