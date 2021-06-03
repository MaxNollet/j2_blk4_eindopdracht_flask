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

    const form_element = document.getElementById("input_form");

    const target_url = form_element.getAttribute("action");
    const request = new XMLHttpRequest();
    request.open("POST", target_url, true);
    // request.setRequestHeader("Content-Type", "application/x-www-form-urlencoded");
    request.onreadystatechange = function () {
        if (this.readyState === 4) {
            let response_message = "";
            let message_type = "";
            switch (this.status) {
                case 200:
                    const results = JSON.parse(this.responseText);
                    response_message = results.message;
                    message_type = results.type;
                    break;
                case 0:
                    response_message = `Could not contact the server! Is the server running and are you ` +
                        `connected to the internet? Please contact your system administrator regarding ` +
                        `this issue with the following status code: <i>HTML status code ${this.status}</i>.`;
                    message_type = "danger"
                    break;
                default:
                    response_message = `Something went wrong on the server! Please contact your system administrator ` +
                        `regarding this issue with the following status code: <i>HTML status code ${this.status}</i>.`
                    message_type = "danger"
                    break
            }
            MessageBuilder("alert_box", message_type, response_message);
            ToggleDisableSubmitButton(false);
        }
    };
    request.send(new FormData(form_element));
}
