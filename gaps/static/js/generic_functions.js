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
