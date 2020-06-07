from Bio import SeqIO
from Bio import pairwise2

fasta_path = '2ABJ.fasta'
seqs = list(SeqIO.parse(fasta_path,'fasta'))
aligns = pairwise2.align.globalxx(seqs[0].seq,seqs[2].seq)
print(pairwise2.format_alignment(*aligns[0]))