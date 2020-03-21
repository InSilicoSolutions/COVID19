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
m.render();