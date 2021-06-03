/**
 * A function which changes the target-attribute of
 * the input form. Switches between 'default' and
 * '_blank'.
 *
 * @author Max Nollet
 * */
function OpenNewTab(source, target) {
    const element_form = document.getElementById(target);
    if (document.getElementById(source).checked) {
        element_form.setAttribute("target", "_blank");
    } else {
        element_form.removeAttribute("target");
    }
}

/**
 * A function which clears the selected file from
 * the file-chooser in the interface if a file has
 * been chosen by the user.
 *
 * @author Max Nollet
 * */
function ClearFile(target) {
    const file_chooser = document.getElementById(target);
    if (file_chooser.value.trim() !== "") {
        file_chooser.value = "";
    }
}

/**
 * A function which builds an alert with a given message and
 * desired type and puts the alert in the specified target
 * element.
 *
 * @author Max Nollet
 * @param target Target element to put the alert into.
 * @param type Type of the alert (success, warning or danger).
 * @param message Message to display in the alert.
 * */
function MessageBuilder(target, type, message) {
    let alert_header = "";
    let type_safe = "";
    switch (type) {
        default:
        case "info":
            alert_header = "Info";
            type_safe = "info";
            break;
        case "success":
            alert_header = "Success";
            type_safe = "success";
            break;
        case "warning":
            alert_header = "Alert";
            type_safe = "warning";
            break;
        case "danger":
            alert_header = "Warning";
            type_safe = "danger";
            break;
    }
    const alert_box = document.getElementById(target);
    const alert_element = document.createElement("div");
    alert_element.setAttribute("class", `alert alert-${type_safe} alert-dismissible`);
    alert_element.setAttribute("role", "alert");
    alert_element.setAttribute("id", `alert_${type_safe}`);
    const alert_dismiss_element = document.createElement("button");
    alert_dismiss_element.setAttribute("type", "button");
    alert_dismiss_element.setAttribute("class", "btn-close");
    alert_dismiss_element.setAttribute("data-bs-dismiss", "alert");
    alert_dismiss_element.setAttribute("aria-label", "Close");
    alert_element.appendChild(alert_dismiss_element);
    const alert_message_element = document.createElement("span");
    alert_message_element.setAttribute("id", `alert_${type_safe}_message`);
    alert_message_element.innerHTML = `<strong>${alert_header}</strong> ${message}`;
    alert_element.appendChild(alert_message_element);
    alert_box.appendChild(alert_element);
}

/**
 * A function which removes all child elements from a
 * given target.
 *
 * @author Max Nollet
 * @param target Target element ot remove all childs from.
 * */
function RemoveChilds(target) {
    const target_element = document.getElementById(target);
    while (target_element.firstChild) {
        target_element.removeChild(target_element.firstChild)
    }
}
