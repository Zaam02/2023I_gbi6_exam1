import csv
from Bio import Entrez
from Bio import SeqIO
from collections import Counter

#------------------2.1---------------------------------------------
def source(email, search_term, output_file):
    Entrez.email = email
    handle = Entrez.esearch(db="nucleotide", term=search_term)
    record = Entrez.read(handle)
    id_list = record['IdList']
    handle = Entrez.efetch(db="nucleotide", id=id_list, rettype="fasta", retmode="text")
    sequences = SeqIO.parse(handle, "fasta")

    organism_names = []
    for sequence in sequences:
        description_parts = sequence.description.split()
        organism_name = " ".join(description_parts[1:3])
        organism_names.append(organism_name)
    species_counts = Counter(organism_names)

    species_counts = Counter(organism_names)

    with open(output_file, "w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(["Species", "Count"])
        for species, count in species_counts.items():
            writer.writerow([species, count])
    
email = "camila.zamora@est.ikiam.edu.ec"
search_term = "Shigel   la"
output_file = "results/source.csv"

source(email, search_term, output_file)

