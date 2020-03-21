import os
import platform
import subprocess
from Bio import GenBank
import time
import dateutil.parser
import requests
import xml.dom.minidom as minidom
import time

validList = ['ORF7B', 'ORF1A', 'ORF10', 'NUCLEOCAPSID', 'ORF8', 'ORF7A', 'ORF6', 'MEMBRANE', 'ENVELOPE', 'ORF3A', 'SURFACE', 'ORF1AB']

geneAlias = {'ORF3':'ORF3A', 'ORF7':'ORF7A', 'N':'NUCLEOCAPSID', 'S':'SURFACE', 'SPIKE':'SURFACE', 'M':'MEMBRANE', 'E': 'ENVELOPE'}

#load reference sequence of each virus gene.  Return it in a hash keyed by gene.
#Put the sequence in a list becasue we will add more sequences when samples are loaded.
def load_ref():
    ref_file = open("data/refernece_sequence.fasta", "r")
    line = ref_file.readline()
    gene = ""
    seq_dict = {}
    while (line):
        if line.startswith(">"):
            if gene != "" and gene in validList:
                seq_dict[gene] = ['Reference|' + sequence]
            gene = line.split()[1].upper()
            sequence = ""
        else:
            sequence += line.strip()
    
        line = ref_file.readline()
        
    if gene in validList:    
        seq_dict[gene] = ['Reference|' + sequence]
    return(seq_dict)    

def findFeature(record, key):
    for feature in record.features:
        if (feature.key == key):
            return feature
        
    return None   

def findAllFeature(record, key):
    fList = []
    for feature in record.features:
        if (feature.key == key):
            fList.append(feature)
        
    return fList   

def findItem(feature, item_name):
    for qualifier in feature.qualifiers:
        if qualifier.key == item_name:
            return qualifier.value.replace('\"','').strip()
        
    return None     

#Load the sequencing samples in the downloaded genbank records.  For now,
#just look for translated protein sequences that identify the gene.  Don't
#worry about the DNA sequences.  Some records have the full virus sequence
#with all the genes.  Some are partial.  Not all records have the /gene 
#annotation with the translation.  For now if they don't have the gene 
#with protein translation I am skipping them.  I am also not trying to handle
#splitting up the big poly protein gene.
def load_samples(sequences):

    with open("data/genbank_sequences.gb") as handle:
        #Use biopython to parse the GenBank records
        for record in GenBank.parse(handle):
            #skip partial sequence records
            if ('partial' in record.definition):
                continue
            
            #For now id for a sample will include a truncated version of the country of origin, date collected
            #and the accession number of the record.
            accession = record.accession[0]
            source = findFeature(record, 'source')
            if source is not None:
                country = findItem(source, '/country=')
                col_date = findItem(source, '/collection_date=')
            id = accession
            if col_date is not None:
                dt = dateutil.parser.parse(col_date) # Time formatting is not consistent
                norm_date = dt.strftime(r'%Y-%m-%d')
                id = norm_date + '-' + id
            if country is not None:
                country = country.replace(':', ' ')
                id = country.split()[0][:7].strip() + '-' + id   
        
            #For each CDS record
            genes = findAllFeature(record, 'CDS')
            for gene in genes:
                #First figure out the gene / protein name.  Try both the product tag and gene tag
                #NOTE: We are not interested in the post translation non structural protein products.
                #      We just process the orf1ab gene that has all of them embedded in it.
                product = findItem(gene, '/product=')
                if product is not None:
                    gene_name = product.split()[0].upper()
                    #A few proteins have aliases, map them to the standard form.
                    if gene_name in geneAlias:
                        gene_name = geneAlias[gene_name]
                
                if gene_name is None or gene_name not in validList:
                    gene_name = findItem(gene, '/gene=') 
                    if (gene_name is not None and gene_name in geneAlias):
                        gene_name = geneAlias[gene_name]
                       
                if (gene_name not in validList ):        
                    continue
                
                sequence = findItem(gene, '/translation=')  
                if (id is not None and sequence is not None):
                    sequences[gene_name].append(id+'|' + sequence)
 

#Create a webpage for a protein to show its multi sequence alignment
#embed the alignment into the page as a string to avoid cross domain errors.
def write_webpage(protein, clustal_file, web_file):
    with open('embed.js') as f:
        embed = f.read()
    out = open(web_file,'w')
    clustal = open(clustal_file, 'r')
    out.write('<meta charset="UTF-8">\n')
    out.write('<html><head>\n')
    out.write('<script src="msa.min.gz.js"></script>\n')
    out.write('</head>\n')
    out.write('<body>\n')
    out.write('<div id="menuDiv"></div>\n')
    out.write('<div id="yourDiv">Loading ... </div>\n')
    out.write('<script>\n')
    out.write('var clustal = \'')
    lines = clustal.readlines()
    for line in lines:
        out.write(line.rstrip() +'\\n')
    out.write('\'\n')
    out.write(embed)
    if not(embed.endswith('\n')):
        out.write('\n')
    out.write('</script>\n')
    out.write('</body></html>\n')
    out.close()
         
def write_index(seq_names, index_fp):
    seq_names.sort()
    with open(index_fp,'w') as out:
        out.write('<meta charset="UTF-8">\n')
        out.write('<html><head><title>COVID19 MSA</title></head>\n')
        out.write('<body>\n')
        out.write('<h1>COVID-19 Sequence Alignments</h1>\n')
        out.write(f'<h4>Last aligned: {time.ctime()}</h4>\n')
        out.write('<ul>\n')
        for seq in seq_names:
            fn = seq+'.html'
            out.write(f'<li><a href="{fn}" target="_blank">{seq}</a></li>\n')
        out.write('</ul>\n')
        out.write('</body>\n')
        out.write('</html>\n')
    
def get_genbank ():
    r = requests.get('http://www.ncbi.nlm.nih.gov/genbank/sars-cov-2-seqs/#nucleotide-sequences')
    dom1 = minidom.parseString(r.text)
    for a in dom1.getElementsByTagName('a'):
        if a.getAttribute('title') == 'SARS-CoV-2 nucleotide sequences':
            link = a
            break
    href = link.getAttribute('href')
    uids = href[href.index('nuccore') + 8:].split(',')
    fn = os.path.join('data', 'genbank_sequences.gb')
    wf = open(fn, 'w')
    print(f'Fetching {len(uids)} GenBank records...')
    for i in range(len(uids)):
        uid = uids[i]
        print(f'    {i + 1}: {uid}')
        r2 = requests.get(f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={uid}&&rettype=gb&retmode=text')
        wf.write(r2.text)
        time.sleep(1)
    wf.close()

get_genbank()
sequences = load_ref()
load_samples(sequences)

#write out one fast file per virus gene with reference sequence and then
#all the collected sample sequences.
for key in sequences:
    fasta_file = "data/"+key+".fasta"
    out = open(fasta_file,'w')
    seqs = sequences[key]
    for seq in seqs:
        toks = seq.split('|')
        out.write('>' + toks[0] + '\n')
        out.write(toks[1] + '\n')
    out.close()
    
    #run clustalo to do a multisequence alignment of each gene
    clustal_output = "data/"+key+".clustal"
    cmd = ['.' + os.sep + 'clustal' + os.sep + 'clustalo', '-i', fasta_file, '-o', clustal_output, '--outfmt=clu', '--force']
    print(' '.join(cmd))
    pl = platform.platform()
    if pl.startswith('Windows'):
        subprocess.run(cmd, shell=True)
    else:
       subprocess.run(cmd)
        
    #write a webpage that dispalys the alignment
    web_page = "data/"+key+".html"
    write_webpage(key, clustal_output, web_page)    
    write_index(list(sequences.keys()), 'data/index.html')
