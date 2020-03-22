from collections import defaultdict
from Bio import SeqIO
import numpy as np
# dict_genome = SeqIO.to_dict(SeqIO.parse("rmark3.fa", "fasta"))

genome = "ABCABCAABABCCACACABCACABA"
seq = ["AB", "BC"]


def kmp_matcher(T , P):
  ''' 
    kmp string matcher algorithme O(len(T) +  len(P))
  '''
  pi, ret, j = prefixfunction(P), [], 0
    
  for i in range(len(T)):
      while j > 0 and T[i] != P[j]:
          j = pi[j - 1]
      if T[i] == P[j]: j += 1
      if j == len(P): 
          ret.append(i - (j - 1))
          j = pi[j - 1]
      
  return ret



def prefixfunction(pattern):
  ''' this compute the pi function in the knp algorithm. pi is an array that contain the neccesary strid to do at any level of the comparaison in the pattern P
  O(len(pattern))
  
  '''
  pi = [0]
        
  for i in range(1, len(pattern)):
      j = pi[i - 1]
      while j > 0 and pattern[j] != pattern[i]:
          j = pi[j - 1]
      pi.append(j + 1 if pattern[j] == pattern[i] else j)
  return pi

print( kmp_matcher(genome , 'ABA'))
#print(scan_genome_given_seq(genome, seq))


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

  def scan_genome_given_seq(self, genome, motif_list, window_lenght=200, stride=1, thrsehold=2):

    """
      input :
        - genome : string representing the long string which is the genome
        - motif_list : list of "motif" to check if it exists in a window
        - window_lenght : length of the window where we must check the "motifs"
        - stride : window displacement stride

      Returns the window whose number of sequences present is greater than a certain threshold.
    """

    list_result = list()

    lenght_result = ((len(genome) - window_lenght) / stride) + 1

    for i in range(int(lenght_result)):

      nb_motif_checked = 0

      for motif in motif_list:
        if motif in genome[i:window_lenght + i]:
          nb_motif_checked += 1

      if thrsehold == nb_motif_checked:
        list_result.append((i, window_lenght + i))

    return list_result




if __name__=='__main__':

  g = Genome('rmark3.fa')
  g.read_fasta()
  genome = g.dict_genome['rmark1'].seq

  #print(genome)
  dcount = defaultdict(int)

  seq = ["GAG", "AUG"]
  for pattern in seq:
    dcount[pattern] = len(kmp_matcher(genome , pattern) )
  
  print(dcount)