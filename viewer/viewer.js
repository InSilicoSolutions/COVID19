var proteinWidget = {};
proteinWidget.ProteinViewer = class {
    constructor (div) {
        this.div = div;
        this.viewer = $3Dmol.createViewer(this.div.id);
        this.residues = {};
        // this.baseStyle = {cartoon:{colorscheme:'chainHetatm'}};
        this.baseStyle = {cartoon:{}};
        this.highlightStyle = {cartoon:{color:'red'}};
        this.alignments = {};
    }
    
    async fetchPDB (pdbId) {
        const response = await fetch(`https://files.rcsb.org/download/${pdbId}.pdb`);
        if (response.ok) {
            return await response.text();
        } else {
            throw(new Error(response.statusText));
        }
    }
    
    async loadPDB (pdbId) {
        const pdbText = await this.fetchPDB(pdbId);
        this.loadPDBText(pdbText);
    }
    
    loadPDBText (pdbText) {
        this.viewer.clear()
            .addModel(pdbText,'pdb')
        this.viewer.setStyle(this.baseStyle)
            .zoomTo()
            .render();
        this._parseModel();
    }

    _parseModel () {
        const atoms = this.viewer.selectedAtoms({})
        // Get unique chains
        const chains = atoms
            .map(atom=>atom.chain)
            .filter((chain,i,all)=>all.indexOf(chain)===i);
        for (let i=0; i<chains.length; i++) {
            // Get unique resis
            const chain = chains[i];
            this.residues[chain] = atoms
                .filter(atom=>atom.chain===chain)
                .map(atom=>atom.resi)
                .filter((resi,i,all)=>all.indexOf(resi)===i);
        }
    }

    highlightChain (chain, start, end) {
        this.viewer.setStyle(
            {
                chain: chain,
                resi: this.residues[chain].slice(start,end)
            },
            this.highlightStyle
        ).render();
    }

    setBaseStyle () {
        this.viewer.setStyle(this.baseStyle).render();
    }

    addAlignment (chain, alignText) {
        let alignId=Math.random().toString(36).substring(7);
        while (this.alignments[alignId]!==undefined){
            alignId=Math.random().toString(36).substring(7);
        }
        this.alignments[alignId] = {
            chain: chain,
            alignment: new proteinWidget.Alignment(alignText)
        }
        return alignId;
    }

    highlight (alignId, beg, end) {
        const chain = this.alignments[alignId].chain;
        const alignment = this.alignments[alignId].alignment;
        const slices = alignment.mapSlice(beg,end);
        console.log(slices);
        for (let i=0; i<slices.length; i++) {
            this.highlightChain(chain, slices[i][0], slices[i][1]);
        }
    }
}

proteinWidget.Alignment = class {
    constructor (alignText) {
        this.alignText = alignText;
        this.mask = []
        const [seqBase, seqSymb, seqTarg] = this.alignText.split('\n').slice(0,3);
        let inAlign = false;
        let indexBase = 0;
        let indexTarg = 0;
        let startBase;
        let startTarg;
        let endBase;
        let endTarg;
        for (let i=0; i<seqSymb.length; i++) {
            if (seqSymb[i] === '|') {
                if (!inAlign) {
                    inAlign=true
                    startBase=indexBase;
                    startTarg=indexTarg;
                    endBase=null;
                    endTarg=null;
                }
                indexBase += 1
                indexTarg += 1
            } else {
                if (inAlign) {
                    inAlign=false;
                    endBase = indexBase;
                    endTarg = indexTarg;
                    this.mask.push([[startBase,endBase],[startTarg,endTarg]]);
                    startBase=null;
                    startTarg=null;
                    endBase=null;
                    endTarg=null;
                }
                indexBase += seqBase[i].match(/[a-z]/i) ? 1 : 0;
                indexTarg += seqTarg[i].match(/[a-z]/i) ? 1 : 0;
            }
        }
        if (inAlign) {
            endBase = indexBase;
            endTarg = indexTarg;
            this.mask.push([[startBase,endBase],[startTarg,endTarg]]);
        }
    }
    mapSlice = (beg, end) => {
        let targRanges = [];
        for (let segIndex=0; segIndex<this.mask.length; segIndex++) {
            let [base, targ] = this.mask[segIndex];
            if (end <= base[0]) { // segment is fully beyond slice
                break;
            } else if (base[1] <= beg) { // segment is fully before slice
                continue;
            } else { // Some amount of overlap
                let normBeg = base[0]<=beg ? beg-base[0] : 0;
                let normEnd = end-base[0];
                targRanges.push([
                    targ[0]+normBeg,
                    Math.min(targ[0]+normEnd, targ[1])
                ]);
            }
        }
        return targRanges;
    }
}