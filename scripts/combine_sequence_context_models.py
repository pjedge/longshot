from __future__ import print_function
# Author : Peter Edge
# Email  : pedge@eng.ucsd.edu

# imports
from collections import defaultdict
import sys

def combine_models(modelfiles):

    count_dict = defaultdict(int)

    for modelfile in modelfiles:
        with open(modelfile,'r') as mf:
            for line in mf:
                el = line.strip().split()
                assert(len(el) == 4)
                key = tuple(el[:3])
                value = int(el[3])
                count_dict[key] += value

    for k,v in count_dict.items():
        line = "{} {} {} {}".format(*k,v)
        print(line)

if __name__ == '__main__':

    if len(sys.argv) < 3 or '-h' in sys.argv or '--help' in sys.argv:
        print('''
        Combines the data for multiple sequence context models created using
        --sequence_context_model option in Longshot into a single model.
        This could be used, for example, to combine models created using different
        chromosomes.

        usage: python combine_sequence_context_models.py [model files] > combined.txt
        ''')
        sys.exit(1)

    combine_models(sys.argv[1:])
