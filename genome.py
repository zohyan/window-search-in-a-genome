class Genome:

  def __init__(self, file):

    self.file = file

  def read_fasta(self, fp):

    name, seq = None, []

    for line in fp:

      line = line.rstrip()

      if line.startswith(">"):
        if name: yield (name, ''.join(seq))
        name, seq = line, []
      else:
        seq.append(line)

    if name:
      yield (name, ''.join(seq))

  def fasta_parse(self, fp):

    key = ''
    for line in fp:
      if line.startswith('>'):
        if key:
          yield key, val
        key, val = line[1:].rstrip().split()[0], ''
      elif key:
        val += line.rstrip()
    if key:
      yield key, val

  def read_genome(self):

    '''
    with open(self.file) as fp:
      for name, seq in self.read_fasta(fp):
        print(seq)
    '''

    print('\n'.join('%s: %s' % keyval for keyval in self.fasta_parse(self.file)))

  def scan_genome(self):

    with open(self.file) as fp:
      for name, seq in self.read_fasta(fp):
        print(len(seq))


if __name__ == "__main__":
  g = Genome('rmark3.fa')
  g.read_genome()


