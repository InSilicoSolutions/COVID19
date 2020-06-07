var v;
const align = `
VVGTFKAKDLIVTPATILKEKPDPNNLVFGTVFTDHMLTVEWSSEFGWEKPHIKPLQNLSLHPGSSALHYAVELFEGLKAFRGVDNKIRLFQPNLNMDRMYRSAVRATLPVFDKEELLECIQQLVKLDQEWVPYSTSASLYIRPAFIGTEPSLGVKKPTKALLFVLLSPVGPYFSSGT------FNPVSLWANPKYVRAWKGGTGDCKMGGNYGSSLFAQCEDVDNGCQQVLWLYGRDHQITEVGTMNLFLYWINEDGEEELATPPLDGIILPGVTRRCILDLAHQWGEFKVSERYLTMDDLTTALEGNRVREMFSSGTACVVCPVSDILYKGETIHIPTMENGPKLASRILSKLTDIQYGREESDWTIVLS
  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||            ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
--GTFKAKDLIVTPATILKEKPDPNNLVFGTVFTDHMLTVEWSSEFGWEKPHIKPLQNLSLHPGSSALHYAVELFEGLKAFRGVDNKIRLFQPNLNMDRMYRSAVRATLPVFDKEELLECIQQLVKLDQEWVPYSTSASLYIRPAFIGTEPSLGVKKPTKALLFVLLSPVGP------XXXXXXFNPVSLWANPKYVRAWKGGTGDCKMGGNYGSSLFAQCEDVDNGCQQVLWLYGRDHQITEVGTMNLFLYWINEDGEEELATPPLDGIILPGVTRRCILDLAHQWGEFKVSERYLTMDDLTTALEGNRVREMFSSGTACVVCPVSDILYKGETIHIPTMENGPKLASRILSKLTDIQYGREESDWTIVLS
  Score=358
`.trim();
var alignId;

window.onload = () => {
  const div = document.querySelector('#viewer');
  v = new proteinWidget.ProteinViewer(div);
  v.loadPDB('2ABJ').then(()=>{
    alignId = v.addAlignment('G',align);
    v.highlight(alignId,165,200);
  });
  
}

const tour = () => {
  const sleep = ms => new Promise( r => setTimeout(r, ms));
  (async () => {
    console.log('Load protein');
    await v.loadPDB('2ABJ')
    await sleep(5000);
    console.log('Add alignment and highlight segment');
    alignId = v.addAlignment('G',align);
    v.highlight(alignId,165,200);
    await sleep(5000);
    console.log('Draw surface');
    v.toggleSurface();
    await sleep(5000);
    console.log('Remove highlight');
    v.unHighlight();
    await sleep(5000);
    console.log('Remove surface');
    v.toggleSurface();
    await sleep(5000);
    console.log('Switch to sphere');
    v.presetStyle('sphere');
    await sleep(5000);
    console.log('Draw surface again');
    v.toggleSurface();
    await sleep(5000);
    console.log('Highlight a different portion');
    v.highlight(alignId,10,200);

  })()

}
