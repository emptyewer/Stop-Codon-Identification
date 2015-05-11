import re, sys, time
import types
from Bio import SeqIO
from Bio import Entrez
from Bio import UniGene
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import httplib
from xenopus.models import Gene

Entrez.email = "vkrishnamani@healthcare.uiowa.edu"

#Get protein record from xenopus leavis or protein with largest similarity score
def get_protsim(record):
    xenopus_protein = {}
    largest_similarity = 0.0
    for prot in record.protsim:
        prot = str(prot)
        protsim_dict = {}
        elements = prot.split(";")
        key = ''
        for element in elements:
            element = element.lstrip()
            split = element.split("=")
            protsim_dict[split[0]] = split[1]
            if split[0] == 'ORG':
                key = split[1]
            if split[0] == 'PCT':
                if largest_similarity < float(split[1]):
                    largest_similarity = float(split[1])
        xenopus_protein[key] = protsim_dict

    if '8355' in xenopus_protein.keys():
        return xenopus_protein['8355']
    else:
        for key in xenopus_protein.keys():
            if float(xenopus_protein[key]['PCT']) == largest_similarity:
                return xenopus_protein[key]

records = SeqIO.parse("XL_cDNA_5660.txt", "fasta")
match_unigene_id = re.compile("/ug=(.+)")
est_records_count = 0
unigene_ids = []
for record in records:
    description = record.description.split()
    for text in description:
        m = match_unigene_id.match(text)
        if m:
            unigene_ids.append(m.group(1))
            est_records_count += 1

unique_unigene_ids = list(set(unigene_ids))
# print "ESTs Count: ", len(unigene_ids)
# print "Unique Count: ", len(unique_unigene_ids)

def get_protein_record(protein_id):
    handle = Entrez.efetch(db="protein", id=protein_id, rettype="gb", retmode="xml")
    protein_record = Entrez.read(handle)
    handle.close()
    return protein_record

def get_nucleotide_record(gene_id):
    handle = Entrez.efetch(db="nucleotide", id=gene_id, rettype="gb", retmode="xml")
    gene_sequence_record = Entrez.read(handle)
    handle.close()
    return gene_sequence_record

input = open("Xl_cleaned.data")
records = UniGene.parse(input)
for record in records:
    if record.ID in unique_unigene_ids:
        db_gene = Gene()
        print "\n", "="*30
        print "Processing UniGene ID : ", record.ID, "\n"
        db_gene.unigene_id = record.ID
        protsim = get_protsim(record)
        if type(protsim) != types.NoneType:
            protein_record = ''
            try:
                protein_record = get_protein_record(protsim['PROTID'])
            except httplib.BadStatusLine:
                protein_record = get_protein_record(protsim['PROTID'])
            print "Protein Similarity Match : ", protsim['PCT'], "%"
            db_gene.protein_similarity = protsim['PCT']
            print "Protein ID : ", protsim['PROTID']
            db_gene.protein_id = protsim['PROTID']
            print "Protein Length : ", protein_record[0]['GBSeq_length']
            db_gene.protein_length = protein_record[0]['GBSeq_length']
            print "Protein Source Organism : ", protein_record[0]['GBSeq_source']
            db_gene.organism = protein_record[0]['GBSeq_source']
            print "Protein Name : ", protein_record[0]['GBSeq_definition']
            db_gene.protein_name = protein_record[0]['GBSeq_definition']
            db_gene.protein_sequence = protein_record[0]['GBSeq_sequence'].upper()
            print " "
            print "Gene ID : ", protein_record[0]['GBSeq_source-db'].split()[-1]
            db_gene.gene_id = protein_record[0]['GBSeq_source-db'].split()[-1]
            gene_sequence_record = ''
            try:
                gene_sequence_record = get_nucleotide_record(protein_record[0]['GBSeq_source-db'].split()[-1])
            except httplib.BadStatusLine:
                gene_sequence_record = get_nucleotide_record(protein_record[0]['GBSeq_source-db'].split()[-1])
            print "Gene Name : ", gene_sequence_record[0]['GBSeq_definition']
            db_gene.gene_name = gene_sequence_record[0]['GBSeq_definition']
            print "Gene Type : ", gene_sequence_record[0]['GBSeq_moltype']
            db_gene.gene_type = gene_sequence_record[0]['GBSeq_moltype']
            print "Gene Length : ", gene_sequence_record[0]['GBSeq_length']
            db_gene.gene_length = gene_sequence_record[0]['GBSeq_length']
            # print "Gene Topology :", gene_sequence_record[0]['GBSeq_topology']
            db_gene.gene_sequence = gene_sequence_record[0]['GBSeq_sequence']
            gene_sequence = gene_sequence_record[0]['GBSeq_sequence']
            for level1 in gene_sequence_record[0]['GBSeq_feature-table']:
                for key in level1.keys():
                    if key == 'GBFeature_key':
                        if level1['GBFeature_key'] == 'CDS':
                            start = int(level1['GBFeature_intervals'][0]['GBInterval_from'])
                            end = int(level1['GBFeature_intervals'][0]['GBInterval_to'])
                            diff = end - start + 1
                            print "Start Codon : ", start
                            db_gene.gene_start_nucleotide = start
                            print "End Codon : ", end
                            db_gene.gene_end_nucleotide = end
                            print "Coding Sequence Length : ", diff
                            db_gene.gene_length = diff
                            start_codon = gene_sequence[start-1:start+2]
                            start_aa = Seq(start_codon, generic_dna).translate()
                            db_gene.gene_start_codon = start_codon
                            db_gene.protein_start_aa = start_aa
                            stop_codon = gene_sequence[end-3:end]
                            stop_aa = Seq(stop_codon, generic_dna).translate()
                            db_gene.gene_stop_codon = stop_codon
                            db_gene.protein_stop_aa = stop_aa
                            print "Start Codon : ", start_codon, "(", start_aa , ")"
                            print "Stop Codon : ", stop_codon, "(", stop_aa , ")"
        else:
            db_gene.gene_sequence = "TRANSCRIBED LOCUS BUT NO MATCHING GENE FOUND"
            db_gene.protein_sequence = "TRANSCRIBED LOCUS BUT NO MATCHING GENE FOUND"
            print "TRANSCRIBED LOCUS BUT NO MATCHING GENE FOUND"
        db_gene.save()

