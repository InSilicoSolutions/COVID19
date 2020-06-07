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