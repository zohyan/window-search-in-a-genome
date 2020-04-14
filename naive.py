from genome import Genome
from time import time

def scan_genome_given_seq(genome, seq, window_lenght=3, stride=1, thersehold=2):

    """
      input :
        - genome : string representing the long string which is the genome
        - motif_list : list of "motif" to check if it exists in a window
        - window_lenght : length of the window where we must check the "motifs"
        - stride : window displacement stride

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

if __name__ == "__main__":
    g = Genome('rmark3.fa')
    g.read_fasta()
    genome = g.dict_genome['rmark1'].seq
    sequence = ["AUG"]

    t1 = time()

    for pattern in sequence:
        print(scan_genome_given_seq(genome, pattern))

    t2 = time()

    print('Elapsed time is %f seconds.' % (t2 - t1))

