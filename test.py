from Bio import SeqIO

# dict_genome = SeqIO.to_dict(SeqIO.parse("rmark3.fa", "fasta"))

sequence = "EWAEWAEW"
substring = ["EWA"]

def find_indexes(seq, sub, window=3):

    result = dict()

    for s in sub:

        intermediate_list = []

        for i in range(0, len(seq), window):
            if seq[i:len(s) + i] == s:
                intermediate_list.append((s, i, len(s)+i))

        result[s] = intermediate_list

    return result

print(find_indexes(sequence, substring))
