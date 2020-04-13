from genome import Genome


def kmp_matcher(T , P):

    pi, ret, j = prefixfunction(P), [], 0

    for i in range(len(T)):
        while j > 0 and T[i] != P[j]:
            j = pi[j - 1]
        if T[i] == P[j]:
            j += 1
        if j == len(P):
            ret.append((i - (j - 1),  i - (j - 1) + len(P)))
            j = pi[j - 1]

    return ret

def prefixfunction(pattern):

    pi = [0]

    for i in range(1, len(pattern)):
        j = pi[i - 1]
        while j > 0 and pattern[j] != pattern[i]:
            j = pi[j - 1]
        pi.append(j + 1 if pattern[j] == pattern[i] else j)
    return pi


if __name__ == "__main__":
    g = Genome('rmark3.fa')
    g.read_fasta()
    genome = g.dict_genome['rmark1'].seq
    sequence = "AUG"

    print(kmp_matcher(genome, sequence))
