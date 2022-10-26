from Bio import SeqIO, Seq
from scipy.spatial import distance_matrix
import pandas as pd
import glob

'''Globals'''
START = "ATG"
STOPS = ["TAA", "TAG", "TGA"]

# fasta_sequences = SeqIO.parse(open("data/bacterial1.fasta"),'fasta')
# for fasta in fasta_sequences:
#     name, sequence = fasta.id, str(fasta.seq)
#     codons = [sequence[i:i+3] for i in range(0,len(sequence),3)]

#     print(name)
#     print(sequence)

def fasta_reader(path):
    fasta_sequences = SeqIO.parse(open(path),'fasta')
    for fasta in fasta_sequences:
        name, sequence, reverse = fasta.id, str(fasta.seq), str(fasta.reverse_complement().seq)

        return (name, sequence, reverse)

def find_codon_seq(seq):
    codons = []
    codon = ""
    
    while (len(seq) >= 3):
        if seq[:3] == START:
            while (len(seq) >= 3 and seq[:3] not in STOPS):
                codon+=seq[:3]
                seq = seq[3:]
            codon+=seq[:3]
            seq = seq[3:]
            if (len(codon) > 100):
                codons.append(codon)
            codon = ""
        seq = seq[1:]

    return codons

def fixcodon(codon):
    size = len(codon)
    if size > 3:
        return codon.tail()
    else:
        return codon

def process_all_frames(sequence):
    result = find_codon_seq(sequence)
    sequence = sequence[1:]
    result += (find_codon_seq(sequence))
    sequence = sequence[1:]
    result += (find_codon_seq(sequence))
    return result

def process_with_reverse(sequence, reverse):
    return process_all_frames(sequence) + process_all_frames(reverse)

def codon_dictionary_generator():
    dictionary = {}
    letters = ['A', 'C', 'G', 'T']
    for l1 in letters:
            for l2 in letters:
                for l3 in letters:
                    dictionary[l1+l2+l3] = 0
                    
    return dictionary

def dicodon_dictionary_generator():
    dictionary = {}
    letters = ['A', 'C', 'G', 'T']
    for l1 in letters:
            for l2 in letters:
                for l3 in letters:
                    for l4 in letters:
                        for l5 in letters:
                            for l6 in letters:
                                dictionary[l1+l2+l3+l4+l5+l6] = 0
                    
    return dictionary

def codon_assign(seq_list):
    dictionary = codon_dictionary_generator()
    for seq in seq_list:
        codons = [seq[i:i+3] for i in range(0,len(seq),3)]
        for codon in codons:
            dictionary[codon] += 1

    return dictionary

def dicodon_assign(seq_list):
    dictionary = dicodon_dictionary_generator()
    for seq in seq_list:
        dicodons = [seq[i:i+6] for i in range(0,len(seq),6)]
        for dicodon in dicodons:
            dictionary[dicodon] += 1

    return dictionary

def file_process(file_list):
    names = []
    codon_dic_list = []
    dicodon_dic_list = []

    for file in file_list:
        name, sequence, reverse = fasta_reader(file)

        names.append(name)

        codons = process_with_reverse(sequence, reverse)
        dicodons = list(filter(lambda cod: len(cod)%6 == 0, codons))

        codon_dic = codon_assign(codons)
        codon_dic_list.append(codon_dic)

        dicodon_dic = dicodon_assign(dicodons)
        dicodon_dic_list.append(dicodon_dic)


    return (names, codon_dic_list, dicodon_dic_list)

def main():
    """ Main program """
    files = glob.glob("data/*.fasta")

    names, codons, dicodons = file_process(files)
    
    # Codons
    data = [list(e.values()) for e in codons]
    ctys = names
    df = pd.DataFrame(data, columns=list(codons[0].keys()), index=ctys)
    matrix = pd.DataFrame(distance_matrix(df.values, df.values), index=df.index, columns=df.index)
    print(matrix)
    matrix.to_csv('codon_matrix.csv', header=False)

    # Dicodons
    data = [list(e.values()) for e in dicodons]
    ctys = names
    df = pd.DataFrame(data, columns=list(dicodons[0].keys()), index=ctys)
    matrix = pd.DataFrame(distance_matrix(df.values, df.values), index=df.index, columns=df.index)
    print(matrix)
    matrix.to_csv('dicodon_matrix.csv', header=False)

if __name__ == "__main__":
    main()


'''
https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance_matrix.html
'''