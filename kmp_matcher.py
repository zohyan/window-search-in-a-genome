from Bio import SeqIO
import numpy as np
# dict_genome = SeqIO.to_dict(SeqIO.parse("rmark3.fa", "fasta"))

genome = "ABCABCAABABCCACACABCACABA"
seq = ["AB", "BC"]


def kmp_matcher(text , pattern):
  pi= prefixfunciton(pattern)
  print(pi)
  q = 0
  for i in range(len(text)):
    print(i)
    while q>0 and pattern[q+1] != text[i] : 
      q = pi[q]
    if pattern[q+1] == text[i]:
      q+=1
    if q == len(pattern):
      print("found pattern with shift" , i-len(pattern))
      q = pi[q]
    



def prefixfunciton(pattern):
  ''' this compute the pi function in the knp algorithm. pi is an array that contain the neccesary strid to do at any level of the comparaison in the pattern P
  '''
  pi = [0]*len(pattern)  # np.zeros(  len(pattern) , int )
  
  pi[0]=0
  k=0
  for q in range(1 , len(pattern)):
    k= pi[k]
    while k>0 and pattern[k+1] != pattern[q]:
      k= pi[k]

    if pattern[k+1]==pattern[q]:
      k+=1
    pi[q]=k

  return pi


kmp_matcher(genome , 'ABA')
#print(scan_genome_given_seq(genome, seq))
