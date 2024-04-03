import re
from pathlib import Path

#fasta
FASTA = Path("Human/hs_ref_GRCh38_chr21.fa")
with FASTA.open() as fasta:
    seq = fasta.readlines()

seqlist = []
name = []
count_a = -1
count_b = 0
for i in range(len(seq)):
    if ">" in seq[i]:
        name += re.findall(r"ref\|(.+?)\|", seq[i])
        count_a += 1
        if count_a != 0:
            seqlist.append("".join(seq[count_b + 1:i]))
        count_b = i
    
    else:
        seq[i] = seq[i].rstrip("\n")

    if i == len(seq) - 1:
        seqlist.append("".join(seq[count_b + 1:i]))


#gff
GFF = Path("Human/ref_GRCh38_scaffolds.gff3")
with GFF.open() as gff:
    f = gff.readlines()

#スペースで分割
info = []
for i in range(len(f)):
    info.append(f[i].split("\t"))


gene = []
cds = []
cds_count = 0
for i in range(len(info)):  
    if len(info[i]) == 9:
        if info[i][2] == "CDS":
            cds.append("")
            cds_name = re.findall(r"gene=(.+?)[;|\n]", info[i][8])
            cds_parent = re.findall(r"Parent=(.+?);", info[i][8])
            cds[cds_count] = [cds_name[0], cds_parent[0], info[i][0], info[i][3], info[i][4], info[i][6]]
            cds_count += 1
            
counter = 0
gene.append(cds[0])
for i in range(1, len(cds)):
    if cds[i][0] == cds[i-1][0]:
        if cds[i][1] == cds[i-1][1]:
            cds[i] += cds[i-1][3:]
            counter += 1
        else:
            gene.append("")
    else:
        gene.append("")
    
    gene[i - counter] = cds[i]


gene_name = []
for i in range(len(gene)):
    gene_name.append(gene[i][0])


codon = {
    "TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
    "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L", 
    "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
    "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
    "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
    "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "TAT":"Y", "TAC":"Y", "TAA":"", "TAG":"",
    "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "TGT":"C", "TGC":"C", "TGA":"", "TGG":"W",
    "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G"
    }


gene_id = input()
n = []
m = ""
l = []
for i in range(len(gene_name)):
    if gene_name[i] == gene_id:
        n.append(i)
        m = name.index(gene[i][2])
        l.append(len(gene[i]))


gene_cds_all = [""] * len(n)
amino_all = [""] * len(n)
num_list_all = [""] * len(n)
for i in range(len(n)):
    gene_cds_0 = ""
    num_list = []
    for j in range(3, l[i], 3):
        num_list.append(int(gene[n[i]][j]))
        num_list.append(int(gene[n[i]][j + 1]))
        num_list.sort()
    num_list_all[i] = num_list


    for j in range(0, len(num_list), 2):
        gene_cds_0 += seqlist[m][(num_list[j] - 1):num_list[j + 1]]

    if gene[n[i]][5] == "+":
        gene_cds = gene_cds_0

    else:
        converter = {"A":"T", "T":"A", "G":"C", "C":"G"}
        gene_cds_convert = [""] * len(gene_cds_0)
        for j in range(len(gene_cds_0)):
            gene_cds_convert[j] = converter[gene_cds_0[j]]
        
        gene_cds = "".join(gene_cds_convert)
        gene_cds = gene_cds[::-1]
    
    gene_cds_all[i] = gene_cds

    amino = [""] * len(gene_cds)
    

    for j in range(0, len(gene_cds), 3):
        amino[j] = codon[gene_cds[j:j+3]]
        
    amino_all[i] = "".join(amino)


print(amino_all)

