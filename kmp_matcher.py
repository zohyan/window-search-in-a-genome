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

if __name__=='__main__':
  dcount = defaultdict(int)
  for pattern in seq:
    dcount[pattern] = len(kmp_matcher(genome , pattern) )
  print(dcount)