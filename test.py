from Bio import SeqIO

dict_genome = SeqIO.to_dict(SeqIO.parse("rmark3.fa", "fasta"))

sequence = dict_genome['rmark1'].seq
substring = ["GAG", "AUG"]

def find_indexes(seq, sub):

    result = dict()

    for s in sub:

        intermediate_list = []

        for i in range(0, len(seq), 1):
            if seq[i:len(s) + i] == s:
                intermediate_list.append((s, i, len(s)+i))

        result[s] = intermediate_list

    return result

print(find_indexes(sequence, substring)["AUG"])
