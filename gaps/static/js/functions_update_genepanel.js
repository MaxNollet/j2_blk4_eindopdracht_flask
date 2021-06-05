/**
 * A function which adds several event triggers to
 * elements in the interface.
 *
 * @author Max Nollet
 * */
window.onload = function () {
    document.getElementById("button_clear")
        .addEventListener("click", function () {
            ClearFile("input_upload_genepanel");
        })
    document.getElementById("form_update_genepanel")
        .addEventListener("submit", SubmitGenepanel)
}

function SubmitGenepanel(event) {
    event.preventDefault();
    RemoveChilds("alert_box");
    const original_text = document.getElementById("button_update_genepanel").innerText;
    ToggleDisableSubmitButton(true, "button_update_genepanel", original_text);

    const form_element = document.getElementById("form_update_genepanel");
    const target_url = form_element.getAttribute("action");
    const request = new XMLHttpRequest();
    request.open("POST", target_url, true);
    request.onreadystatechange = function() {
        ResponseHandler(this, "alert_box");
        ToggleDisableSubmitButton(false, "button_update_genepanel", original_text);
    }
    request.send(new FormData(form_element));
}