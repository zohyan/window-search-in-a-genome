from Bio import SeqIO

# dict_genome = SeqIO.to_dict(SeqIO.parse("rmark3.fa", "fasta"))

genome = "EWAEWA"
seq = ["EW", "WA"]

def scan_genome_given_seq(genome, seq, window_lenght=3, stride = 1):

    dict_result = dict()

    lenght_result = ((len(genome) - window_lenght)/stride) + 1

    for s in seq:

        intermediate_list = []
        for i in range(int(lenght_result)):
            if s in genome[i:window_lenght + i]:
                intermediate_list.append((genome[i:window_lenght + i], i, window_lenght + i))

        dict_result[s] = intermediate_list

    return dict_result

print(scan_genome_given_seq(genome, seq))
