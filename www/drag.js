function allowDrop(ev) {
    ev.preventDefault();
}
var y= null

function drop(ev) {
    ev.preventDefault();
    var data = ev.dataTransfer.getData("text");
    Shiny.onInputChange("draggedfile", [y,data])
	y= null

}

