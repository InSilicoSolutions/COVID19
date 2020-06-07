var proteinWidget = {};
proteinWidget.ProteinViewer = class {
    constructor (div) {
        this.div = div;
        this.viewer = $3Dmol.createViewer(this.div.id);
        this.residues = {};
        this.baseStyle = {cartoon:{}};
        this.highlightStyle = {cartoon:{color:'red'}};
        this.alignments = {};
        this.showSurface = false;
        this.surfaceType = 'SAS';
        this.baseSurfaceStyle = {opacity:0.85, color:'grey'}
        this.highlightSurfaceStyle = {opacity:0.85, color:'red'}
        this.highlightSelector = {invert:true};
        this.highlighted = false;
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
        this._parseResidues();
    }

    _parseResidues () {
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

    _drawHighlight () {
        this.viewer.setStyle(
            this.highlightSelector,
            this.highlightStyle
        );
    }

    highlight (alignId, beg, end) {
        const chain = this.alignments[alignId].chain;
        const alignment = this.alignments[alignId].alignment;
        const slices = alignment.mapSlice(beg,end);
        let resis = slices
            .map(slice=>this.residues[chain].slice(slice[0],slice[1]))
            .flat();
        this.highlightSelector = {chain:chain,resi:resis};
        this._drawHighlight();
        this.viewer.render();
        this.highlighted = true;
        this.toggleSurface(this.showSurface);
    }

    toggleSurface (showSurface) {
        showSurface = showSurface===undefined ? !this.showSurface : Boolean(showSurface);
        this.viewer.removeAllSurfaces();
        if (showSurface) {
            if (this.highlighted) {
                this.viewer.addSurface(
                    this.surfaceType,
                    this.baseSurfaceStyle,
                    {...this.highlightSelector,invert:true},
                    {}
                );
                this.viewer.addSurface(
                    this.surfaceType,
                    this.highlightSurfaceStyle,
                    this.highlightSelector,
                    {}
                );
            } else {
                this.viewer.addSurface(
                    this.surfaceType,
                    this.baseSurfaceStyle,
                    {},
                    {}
                )
            }
        }
        this.showSurface = showSurface;
    }
    
    render () {
        this.viewer.setStyle({},this.baseStyle);
        this._drawHighlight();
        this.viewer.render();
    }

    unHighlight() {
        this.highlightSelector = {invert:true}
        this.render();
        this.highlighted = false;
        this.toggleSurface(this.showSurface);
    }
}

proteinWidget.Alignment = class {
    constructor (alignText) {
        this.alignText = alignText;
        this.mask = []
        const [seqBase, seqSymb, seqTarg] = this.alignText
            .split('\n')
            .slice(0,3);
        let inAlign = false;
        let indexBase = 0; // Index on the base sequence, ignoring gaps
        let indexTarg = 0; // Index on the target sequence, ignoring gaps
        let startBase;
        let startTarg;
        let endBase;
        let endTarg;
        for (let i=0; i<seqSymb.length; i++) {
            if (seqSymb[i] === '|') { // In an alignment
                if (!inAlign) { // Start of align
                    inAlign=true
                    startBase=indexBase;
                    startTarg=indexTarg;
                    endBase=null;
                    endTarg=null;
                }
                indexBase += 1
                indexTarg += 1
            } else { // In a diff of some kind
                if (inAlign) { // Start of diff
                    inAlign=false;
                    endBase = indexBase;
                    endTarg = indexTarg;
                    this.mask.push([
                        [startBase,endBase],
                        [startTarg,endTarg]
                    ]);
                    startBase=null;
                    startTarg=null;
                    endBase=null;
                    endTarg=null;
                }
                // Increment indices if sequence is not gapped
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
        for(const section of this.mask){
            let [base, targ] = section;
            if (end <= base[0]) { // section is fully beyond slice
                break;
            } else if (base[1] <= beg) { // section is fully before slice
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