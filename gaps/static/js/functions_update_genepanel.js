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
/**
 * A function which sends an AJAJ-request to the server so
 * that a genepanel can be updated. Possible errors are
 * displayed in an appropriate message-box so the user
 * can act accordingly.
 *
 * @author Max Nollet
 * @param event Default submit-event from submitting the form.
 */
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
        if (request.readyState === 4) {
            ResponseHandler(this, "alert_box");
            ToggleDisableSubmitButton(false, "button_update_genepanel", original_text);
        }
    }
    request.send(new FormData(form_element));
}