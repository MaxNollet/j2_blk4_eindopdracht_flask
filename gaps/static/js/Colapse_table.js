var coll = document.getElementsByClassName("test");
var i;

for (let j = 0; j < coll.length; j++) {
    coll[j].addEventListener("click", function () {
        this.classList.toggle("active");
        var genes = this.nextElementSibling;
        if (genes.style.maxHeight) {
            genes.style.maxHeight = null;
        } else {
            genes.style.maxHeight = genes.scrollHeight + "px";
        }
    })
}
