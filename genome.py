from Bio import SeqIO


class Genome:

  def __init__(self, file):

    self.file = file
    self.dict_genome = None

  def read_fasta(self):
    """
    Reads a fasta file using SeqIO
    """
    self.dict_genome = SeqIO.to_dict(SeqIO.parse(self.file, "fasta"))

  def check_position_sequence_in_genome(self, genome, seq):
    """
    Returns a dictionary containing the start
    and end of each sequence inside genome.
    """
    result = dict()

    for s in seq:

      intermediate_list = []

      for i in range(0, len(genome), 1):
        if genome[i:len(s) + i] == s:
          intermediate_list.append((s, i, len(s) + i))

      result[s] = intermediate_list

    return result

