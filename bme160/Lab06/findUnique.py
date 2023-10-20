#!/usr/bin/env python3
# Name: Riya Tiloda (rtiloda)
# Group Members: Akash Pandit, Natalie Cellucci, Hema Prasanna
import sys
from collections import defaultdict

class FastAreader:
    def __init__(self, fname=''):
        self.fname = fname

    def doOpen(self):
        if self.fname == '':
            return sys.stdin
        else:
            return open(self.fname)

    def readFasta(self):
        with self.doOpen() as fileH:
            header = ''
            sequence = ''
            for line in fileH:
                if line.startswith('>'):
                    if header and sequence:
                        yield header, sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else:
                    sequence += ''.join(line.rstrip().split()).upper()
            if header and sequence:
                yield header, sequence


class tRNA:
    '''For  given input tRNA 'sequence' and 'header', class finds unique subsequences, gets essential subsequences from those unique subsequences. The class makes use of a class attribute alltRNAs.'''

    alltRNAs = []   #class attribute that can be accessed by all tRNA objects

    def __init__(self, sequence):
        tRNA.alltRNAs.append(self)

        s = sequence.replace('-', '').replace('.', '').replace('_', '')

        self.powSet = set()   # set to store all substrings of the sequence
        self.altRNAs = set()  # set to store alternative substrings of other tRNA objects
        
        # Generate all possible substrings of the sequence
        for i in range(len(s)):
            for j in range(i + 1, len(s) + 1):
                self.powSet.add(s[i:j])

    def searchEssential(self):
        '''Find unique elements and return only the essential unique elements'''

        self.altRNAs = set().union(*[rna.powSet for rna in tRNA.alltRNAs if rna is not self])
        self.allUniqes = self.powSet.difference(self.altRNAs)

        self.discardStrings = set()
        for unique in self.allUniqes:
            # Determine if the substring or its left or right (adjacent) neighbor occurs in the unique set
            l = unique[:len(unique) - 1]
            r = unique[1:]
            if l in self.allUniqes or r in self.allUniqes:
                self.discardStrings.add(unique)

        self.essentials = self.allUniqes.difference(self.discardStrings)
        return self.essentials


def main(inCL=None, inFile=None):
    '''Reads a FastA file containing tRNA sequences. Proceeds to print and align unique and essential susbsequences and sorts alphabetically.'''
    reader = FastAreader(inFile)

    # Read and process each FASTA entry
    for header, sequence in sorted(reader.readFasta()):
        seq = sequence.replace('-', '').replace('.', '').replace('_', '')
        tRNAinst = tRNA(seq)

        print(header)
        print(sequence)

        essentialSet = tRNAinst.searchEssential()

        # Print the essential substrings, which are sorted on basis of sequence position
        for item in sorted(essentialSet, key=lambda x: sequence.find(x)):
            print('-' * (sequence.find(item) - 1), item)


if __name__ == "__main__":
    main(inFile='bos-tRNA-7.fa')