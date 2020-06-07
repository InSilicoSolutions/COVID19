class ProteinViewer {
    constructor (div) {
        this.div = div;
        this.viewer = $3Dmol.createViewer(this.div.id);
        this.residues = {};
        this.baseStyle = {cartoon:{}};
        this.highlightStyle = {cartoon:{color:'red'}};
    }
    
    async getPDB (pdbId) {
        const response = await fetch(`https://files.rcsb.org/download/${pdbId}.pdb`);
        if (response.ok) {
            return await response.text();
        } else {
            throw(new Error(response.statusText));
        }
    }
    
    async loadPDB (pdbId) {
        const pdbContent = await this.getPDB(pdbId);
        this.loadPDBFile(pdbContent);
    }
    
    loadPDBFile (pdbContent) {
        this.viewer.clear();
        this.viewer.addModel(pdbContent,'pdb');
        this.viewer.setStyle(this.baseStyle)
        this.viewer.center().zoomTo().render();
        this.parseModel();
    }

    parseModel () {
        const atoms = this.viewer.selectedAtoms({})
        const chains = atoms
            .map(atom=>atom.chain)
            .filter((chain,i,all)=>all.indexOf(chain)===i);
        for (let i=0; i<chains.length; i++) {
            const chain = chains[i];
            this.residues[chain] = atoms
                .filter(atom=>atom.chain===chain)
                .map(atom=>atom.resi)
                .filter((resi,i,all)=>all.indexOf(resi)===i);
        }
    }

    highlightChain (chain, start, end) {
        console.log(chain,start,end);
        this.viewer.setStyle(
            {
                chain:chain,
                resi:this.residues[chain].slice(start,end)
            },
            this.highlightStyle
        ).render();
    }

    setBaseStyle () {
        this.viewer.setStyle(this.baseStyle).render();
    }

    static alignedSlices = (aligns, beg, end) => {
        let targRanges = [];
        for (let segIndex=0; segIndex<aligns.length; segIndex++) {
            let [base, targ] = aligns[segIndex];
            if (end <= base[0]) { //tail
                break;
            } else if (base[1] <= beg) { //head
                continue;
            } else {
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