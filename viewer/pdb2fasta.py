from Bio import SeqIO
from Bio.PDB import PDBList

pdbid = '2ABJ'
pdbdir = 'pdb'
pdbl = PDBList(pdb=pdbdir)
struct_path = pdbl.retrieve_pdb_file(pdbid)
fasta_path = f'{pdbid}.fasta'
SeqIO.convert(struct_path,'cif-atom',fasta_path,'fasta')