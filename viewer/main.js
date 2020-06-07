var v;
window.onload = () => {
    const div = document.querySelector('#viewer');
    v = new ProteinViewer(div);
    v.loadPDB('2ABJ');
}