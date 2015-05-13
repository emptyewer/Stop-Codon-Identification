import re, sys, datetime, os
import httplib, textwrap, types
import subprocess
from shlex import split
from Bio import SeqIO
from Bio import Entrez
from Bio import UniGene
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

tax_id = sys.argv[1]
expression_tag_filename = ''
unigene_db_filename = ''
output_filename = ''

#Add/Change this section for adding more organisms
if tax_id == '8355':
    from xenopus.models import Gene
    expression_tag_filename = "XL_cDNA_5660.txt"
    unigene_db_filename = "Xl.data"
    output_filename = "results_xenopus.txt"
    if os.path.exists(output_filename):
        output_filename = output_filename[-4] + str(datetime.time()) + ".txt"
elif tax_id == '9606':
    from hek293.models import Gene
    expression_tag_filename = "HS_cDNA_16145.txt"
    unigene_db_filename = "Hs.data"
    output_filename = "results_hek293.txt"
    if os.path.exists(output_filename):
        output_filename = output_filename[-4] + str(datetime.time()) + ".txt"


Entrez.email = "vkrishnamani@healthcare.uiowa.edu"

def parse_tmprediction(output):
    split_values = output[0].split("\n")
    string = ''
    n_tm_helices = 0
    segments = []
    for value in split_values:
        if re.match(r'^%pred', value):
            string = value
    if string != '':
        (title, prediction) = string.split(":")
        segments = prediction.split(',')
        segments = [x.lstrip().rstrip() for x in segments]
        for segment in segments:
            if re.match(r'^M', segment):
                n_tm_helices += 1
    return(n_tm_helices, segments)


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

    if tax_id in xenopus_protein.keys():
        return xenopus_protein[tax_id]
    else:
        for key in xenopus_protein.keys():
            if float(xenopus_protein[key]['PCT']) == largest_similarity:
                return xenopus_protein[key]

records = SeqIO.parse(expression_tag_filename, "fasta")
print "Est db", expression_tag_filename
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
# output.write("ESTs Count: ", len(unigene_ids)
# output.write("Unique Count: ", len(unique_unigene_ids)

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
output = open(output_filename, 'w')
unigene_input = open(unigene_db_filename)
print "Unigene db", unigene_db_filename
unigene_records = UniGene.parse(unigene_input)
processed = 1
for record in unigene_records:
    if record.ID in unique_unigene_ids and not Gene.objects.filter(unigene_id=record.ID).exists():
        db_gene = Gene()
        db_gene.taxonomy_id = tax_id
        output.write("="*30)
        output.write("\n")
        output.write("UniGene ID : %s\n" % record.ID)
        print "Accessing Unigene Record... ", record.ID, " (", processed, " of ", est_records_count, ")"
        db_gene.unigene_id = record.ID
        protsim = get_protsim(record)
        if type(protsim) != types.NoneType:
            protein_record = ''
            db_gene.protein_similarity_taxonomy_id = protsim['ORG']
            try:
                protein_record = get_protein_record(protsim['PROTID'])
            except httplib.BadStatusLine:
                protein_record = get_protein_record(protsim['PROTID'])
            output.write("Protein Similarity Match : %s %%\n" % protsim['PCT'])
            db_gene.protein_similarity = protsim['PCT']
            output.write("Protein ID : %s\n" % protsim['PROTID'])
            db_gene.protein_id = protsim['PROTID']
            output.write("Protein Length : %s\n" % protein_record[0]['GBSeq_length'])
            db_gene.protein_length = protein_record[0]['GBSeq_length']
            output.write("Protein Source Organism : %s\n" % protein_record[0]['GBSeq_source'])
            db_gene.source_organism = protein_record[0]['GBSeq_source']
            output.write("Protein Name : %s\n" % protein_record[0]['GBSeq_definition'])
            db_gene.protein_name = protein_record[0]['GBSeq_definition']
            protein_sequence = protein_record[0]['GBSeq_sequence'].upper()
            db_gene.protein_sequence = protein_sequence

            commandline = 'decodeanhmm -f /home/venky/Softwares/bin/tmhmm-2.0c/lib/TMHMM2.0.options -modelfile /home/venky/Softwares/bin/tmhmm-2.0c/lib/TMHMM2.0.model'
            command = split(commandline)
            process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            process.stdin.write(">" + protein_record[0]['GBSeq_definition'] + "\n" + protein_sequence)
            (n_tm, segments) = parse_tmprediction(process.communicate())
            if n_tm > 0:
                db_gene.protein_classification = 'Membrane'
                output.write("Protein Classification : Membrane\n")
                db_gene.protein_tm_helices = n_tm
                output.write("Number of Predicted TM helices: %d\n" % (n_tm))
                db_gene.protein_predicted_segments = ", ".join(segments)
                output.write("Predicted Segments : %s\n" % (", ".join(segments)))
            else:
                db_gene.protein_classification = 'Soluble'
                output.write("Protein Classification : Soluble\n")

            for line in textwrap.wrap(protein_sequence, width=100):
                output.write("%s\n" % line)
            output.write("\n")

            output.write("Gene ID : %s\n" % protein_record[0]['GBSeq_source-db'].split()[-1])
            db_gene.gene_id = protein_record[0]['GBSeq_source-db'].split()[-1]
            gene_sequence_record = ''
            try:
                gene_sequence_record = get_nucleotide_record(protein_record[0]['GBSeq_source-db'].split()[-1])
            except httplib.BadStatusLine:
                gene_sequence_record = get_nucleotide_record(protein_record[0]['GBSeq_source-db'].split()[-1])
            output.write("Gene Name : %s\n" % gene_sequence_record[0]['GBSeq_definition'])
            db_gene.gene_name = gene_sequence_record[0]['GBSeq_definition']
            output.write("Gene Type : %s\n" % gene_sequence_record[0]['GBSeq_moltype'])
            db_gene.gene_type = gene_sequence_record[0]['GBSeq_moltype']
            output.write("Gene Length : %s\n" % gene_sequence_record[0]['GBSeq_length'])
            db_gene.gene_length = gene_sequence_record[0]['GBSeq_length']
            # output.write("Gene Topology :", gene_sequence_record[0]['GBSeq_topology']
            db_gene.gene_sequence = gene_sequence_record[0]['GBSeq_sequence']
            gene_sequence = gene_sequence_record[0]['GBSeq_sequence']
            for level1 in gene_sequence_record[0]['GBSeq_feature-table']:
                for key in level1.keys():
                    if key == 'GBFeature_key':
                        if level1['GBFeature_key'] == 'CDS':
                            start = int(level1['GBFeature_intervals'][0]['GBInterval_from'])
                            end = int(level1['GBFeature_intervals'][0]['GBInterval_to'])
                            diff = end - start + 1
                            output.write("CDS Start : %d\n" % start)
                            db_gene.gene_start_nucleotide = start
                            output.write("CDS End : %d\n" % end)
                            db_gene.gene_end_nucleotide = end
                            output.write("Coding Sequence Length : %d\n" % diff)
                            db_gene.gene_length = diff
                            start_codon = gene_sequence[start-1:start+2]
                            start_aa = Seq(start_codon, generic_dna).translate()
                            db_gene.gene_start_codon = start_codon
                            db_gene.protein_start_aa = start_aa
                            stop_codon = gene_sequence[end-3:end]
                            stop_aa = Seq(stop_codon, generic_dna).translate()
                            db_gene.gene_stop_codon = stop_codon
                            db_gene.protein_stop_aa = stop_aa
                            output.write("Start Codon : %s (%s)\n" % (start_codon, start_aa))
                            output.write("Stop Codon : %s (%s)\n" % (stop_codon,  stop_aa))
            for line in textwrap.wrap(gene_sequence, width=100):
                output.write("%s\n" % line)
        else:
            db_gene.gene_sequence = "TRANSCRIBED LOCUS BUT NO MATCHING GENE FOUND"
            db_gene.protein_sequence = "TRANSCRIBED LOCUS BUT NO MATCHING GENE FOUND"
            output.write("TRANSCRIBED LOCUS BUT NO MATCHING GENE FOUND\n")
        db_gene.save()
        processed += 1
unigene_input.close()
output.close()
