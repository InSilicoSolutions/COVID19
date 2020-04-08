var blob = new Blob([clustal], { type: 'text/plain' });
var file = new File([blob], "SURFACE.clustal", {type: "text/plain"});
var opts = {
    el: document.getElementById("yourDiv"), 
    vis: {
        conserv: false, 
        overviewbox: false, 
        seqlogo: true,
    }, 
    conf: { 
        dropImport: true, 
    }, 
    zoomer: { 
        menuFontsize: "12px", 
        autoResize: true, 
        alignmentHeight: 700,
        alignmentWidth: "auto",
        labelNameLength: 180,
    },
};
var m = new msa.msa(opts);
m.u.file.importFile(clustal);
var schemer = m.g.colorscheme;
var refSeq = m.seqs.get(0).attributes.seq;
schemer.addDynScheme('altered', (letter, info) => {
    if (info.pos < refSeq.length && letter !== refSeq[info.pos]) {
        return '#f00';
    } else if (info.y % 2) {
        return '#aaa';
    } else {
        return '#ddd';
    }
})
schemer.set('scheme','altered');
var menuOpts = {el: document.getElementById('div'), msa: m};
var defMenu = new msa.menu.defaultmenu(menuOpts);
m.addView("menu", defMenu);
reSize();
m.render();

function injectOption(){
	var optionName = "Differences"
	var styleDropdown = document.getElementsByClassName("smenu-dropdown-menu")[5];
	firstOption = styleDropdown.childNodes[0]
	if (firstOption.textContent != optionName){
		var newOption = document.createElement("LI"); 
		newOption.style = "line-height: 14px;";
		newOption.textContent = optionName;
		newOption.onclick = styleSwitch
		styleDropdown.insertBefore(newOption, firstOption);
	}
	
}

function styleSwitch(){
	schemer.set('scheme','altered');
}

function reSize(){
	m.g.zoomer.set("alignmentWidth", window.innerWidth-220)
	m.g.zoomer.set("alignmentHeight", window.innerHeight-183)
}
document.onclick = injectOption
window.addEventListener("resize", reSize, false);