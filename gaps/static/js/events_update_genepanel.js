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
}
