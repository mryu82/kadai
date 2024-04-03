#fasta
fasta = open("S_cerevisiae/saccharomyces_cerevisiae_R64-1-1_20110208_sequences.fasta", "r")

seq = fasta.read()
seqlist = seq.split("\n")
seqlist_new = []
name = []
count = -1

for i in range(len(seqlist)):
    if ">" in seqlist[i]:
        name.append(seqlist[i][1:])
        count += 1
        seqlist_new.append("")
    else:
        seqlist_new[count] += seqlist[i]


#gff

import re

gff = open("S_cerevisiae/saccharomyces_cerevisiae_R64-1-1_20110208_annotation.gff", "r")
f = gff.readlines()

# #で始まる行を削除
for i in range(len(f)):
    if f[0][0] == "#":
        f.pop(0)
    else:
        break

#スペースで分割
inf = []
for i in range(len(f)):
    inf.append("")
    inf[i] = f[i].split("\t")


gene = []
cds = []
cds_count = 0
for i in range(len(f)):  
    
    if inf[i][2] == "CDS":
        cds.append("")
        cds_name = re.findall(r"Parent=(.*?);", inf[i][8])
        cds[cds_count] = [cds_name[0], inf[i][0], inf[i][3], inf[i][4], inf[i][6]]
        cds_count += 1


counter = 0
gene.append(cds[0])
for i in range(1, len(cds)):
    if cds[i][0] == cds[i-1][0]:
        cds[i] += cds[i-1][2:]
        counter += 1
    else:
        gene.append("")
    
    gene[i - counter] = cds[i]


#辞書{遺伝子名：配列}
dict = {}

for i in range(len(gene)):
    m = name.index(gene[i][1])
    l = len(gene[i])
    
    gene_cds_0 = ""
    for j in range(1, (l // 3) + 1):
        gene_cds_0 += seqlist_new[m][(int(gene[i][l - 3*j]) - 1):int(gene[i][l - 3*j + 1])]

    gene_cds = ""
    if gene[i][4] == "+":
        gene_cds = gene_cds_0

    else:
        for j in range(len(gene_cds_0)):
            if gene_cds_0[j] == "A":
                gene_cds += "T"
            if gene_cds_0[j] == "T":
                gene_cds += "A"
            if gene_cds_0[j] == "G":
                gene_cds += "C"
            if gene_cds_0[j] == "C":
                gene_cds += "G"
    
    dict[gene[i][0]] = gene_cds


id = input()
print(dict[id])