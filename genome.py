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

  def scan_genome_naive(self, seq, sub, window=1):
    """
    Returns a dictionary containing the start
    and end of each sub sequence find window.
    """
    result = dict()

    for s in sub:

      intermediate_list = []

      for i in range(0, len(seq), window):
        if seq[i:len(s) + i] == s:
          intermediate_list.append((s, i, len(s) + i))

      result[s] = intermediate_list

    return result

if __name__ == "__main__":
  g = Genome('rmark3.fa')
  g.read_fasta()
  sequence = g.dict_genome['rmark1'].seq
  substring = ["GAG", "AUG"]
  print(g.scan_genome_naive(sequence, substring))


