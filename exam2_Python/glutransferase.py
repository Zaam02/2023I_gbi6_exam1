import csv
from Bio import Entrez
from Bio import SeqIO
from collections import Counter
from Bio.SeqUtils import ProtParam
import matplotlib.pyplot as plt

#------- 2.1 -----
def source(email, input_file, output_file):
    Entrez.email = email

    with open(input_file, 'r') as file:
        accession_numbers = [line.strip() for line, _ in zip(file, range(50))] 
    organism_names = []
    for accession_number in accession_numbers:
        handle = Entrez.efetch(db="nucleotide", id=accession_number, rettype="gb", retmode="text")
        record = SeqIO.read(handle, 'gb')
        description_parts = record.description.split()
        organism_name = " ".join(description_parts[1:3])
        organism_names.append(organism_name)

    species_counts = Counter(organism_names)

    with open(output_file, "w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(["Species", "Count"])
        for species, count in species_counts.items():
            writer.writerow([species, count])

email = "camila.zamora@est.ikiam.edu.ec"
input_file = "data/gstm.txt"
output_file = "results/source.csv"

source(email, input_file, output_file)

#-------------- 2.2 -------------
def sequences():
    Entrez.email = 'camila.zamora@est.ikiam.edu.ec'

    accession_numbers = []
    with open('data/gstm.txt', 'r') as f:
        reader = csv.reader(f)
        for i, row in enumerate(reader):
            if i >= 90:  
                break
            accession_numbers.append(row[0].strip())

    for accession_number in accession_numbers:
        handle = Entrez.efetch(db='nucleotide', id=accession_number, rettype='gb', retmode='text')
        record = SeqIO.read(handle, 'gb')
        dna_sequence = record.seq

        # Traducción y separación de péptidos
        protein_sequence = dna_sequence.translate()
        peptides = protein_sequence.split('*')

        # Filtrar péptidos que comienzan con Metionina
        met_peptides = [peptide for peptide in peptides if peptide.startswith('M')]

        # Calcular el peso molecular e índice de inestabilidad para cada péptido
        peptide_data = []
        for peptide in met_peptides:
            analysis = ProtParam.ProteinAnalysis(str(peptide))
            molecular_weight = analysis.molecular_weight()
            instability_index = analysis.instability_index()
            peptide_data.append((peptide, molecular_weight, instability_index))

            # Guardar los resultados en un archivo CSV
            with open('results/glupeptides.csv', 'a') as f:
                f.write(f'{accession_number},{peptide},{molecular_weight},{instability_index}\n')

        # Crear el gráfico de dispersión
        mw_values = [data[1] for data in peptide_data]
        ii_values = [data[2] for data in peptide_data]

        plt.scatter(mw_values, ii_values, color='#20B2AA', s=90, marker='+')
        plt.xlabel('Peso Molecular')
        plt.ylabel('Índice de estabilidad')
        plt.title('Glucopeptidos - Peso molecular vs Índice de estabilidad')
        plt.savefig('results/glupeptides.png')
        plt.close()

sequences()