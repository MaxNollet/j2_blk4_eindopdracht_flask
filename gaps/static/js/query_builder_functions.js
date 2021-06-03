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

function ToggleDisableSubmitButton(disable) {
    const submit_button = document.getElementById("button_submit");
    const spinner = document.getElementById("spinner_submit_button");
    const button_text = document.getElementById("text_submit_button");
    if (disable === true) {
        submit_button.setAttribute("disabled", "true");
        spinner.removeAttribute("hidden");
        button_text.innerText = "Processing...";
    } else {
        button_text.innerText = "Search genes";
        spinner.setAttribute("hidden", "true");
        submit_button.removeAttribute("disabled");
    }
}

/***/
function Submit(event) {
    event.preventDefault();
    RemoveChilds("alert_box");
    ToggleDisableSubmitButton(true);

    const target_url = document.getElementById("input_form").getAttribute("action");
    const request = new XMLHttpRequest();
    request.onreadystatechange = function () {
        if (this.readyState === 4) {
            if (this.status === 200) {
                console.log(this.responseText)
                MessageBuilder("alert_box", "info", this.responseText)
            } else {
                let error_message = `Unknown error! Please contact your system administrator with this code` +
                    `: <i>HTML status code ${this.status}</i>`;
                switch (this.status) {
                    case 0:
                        error_message = "Could not contact the server! If the server running and are you " +
                            "connected to the internet? <i>(HTML status code 0)</i>";
                        break;
                }
                MessageBuilder("alert_box", "danger", error_message);
            }
            ToggleDisableSubmitButton(false);
        }
    };
    request.open("POST", target_url, true);
    request.send();
}
