from Bio import SeqIO

# dict_genome = SeqIO.to_dict(SeqIO.parse("rmark3.fa", "fasta"))

genome = "ABCABC"
seq = ["AB", "BC"]

def scan_genome_given_seq(genome, seq, window_lenght=3, stride=1, thersehold=2):

    """
    input :
    - genome
    - seq
    - window_lenght
    - stride

    Returns the window whose number of sequences present is greater than a certain threshold.
    """

    list_result = list()

    lenght_result = ((len(genome) - window_lenght)/stride) + 1

    for i in range(int(lenght_result)):

        nb_motif_checked = 0

        for s in seq:
            if s in genome[i:window_lenght + i]:
                nb_motif_checked += 1

        if thersehold == nb_motif_checked:
            list_result.append((i, window_lenght + i))

    return list_result

print(scan_genome_given_seq(genome, seq))
