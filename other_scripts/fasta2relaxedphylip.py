import re
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', metavar='input_fasta', dest='i',
                    type=str, required=True)
parser.add_argument('-o', '--output', metavar='output_phylip', dest='o',
                    type=str, required=True)
args = parser.parse_args()

with open(args.i, 'r') as fin:
    sequences = [(m.group(1), ''.join(m.group(2).split()))
    for m in re.finditer(r'(?m)^>([^ \n]+)[^\n]*([^>]*)', fin.read())]
with open(args.o, 'w') as fout:
    fout.write('%d %d\n' % (len(sequences), len(sequences[0][1])))
    for item in sequences:
        fout.write('%-14s %s\n' % item)
