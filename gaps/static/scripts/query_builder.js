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
    alert("TEST!!")
}

function Concatinate() {

}
