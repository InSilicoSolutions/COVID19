class ProteinViewer {
    constructor (div) {
        this.div = div;
        this.viewer = $3Dmol.createViewer(this.div.id);
        this.residues = {}
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
        this.viewer.setStyle({cartoon:{}})
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
}