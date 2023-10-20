#!/usr/bin/env python3
# Name: Riya Tiloda (rtiloda)
# Group Members: Hema Prasanna, Akash Pandit
import sys
class FastAreader :
    ''' 
    Define objects to read FastA files.
    
    instantiation: 
    thisReader = FastAreader ('testTiny.fa')
    usage:
    for head, seq in thisReader.readFasta():
        print (head,seq)
    '''
    def __init__ (self, fname=None):
        '''contructor: saves attribute fname '''
        self.fname = fname
            
    def doOpen (self):
        ''' Handle file opens, allowing STDIN.'''
        if self.fname is None:
            return sys.stdin
        else:
            return open(self.fname)
        
    def readFasta (self):
        ''' Read an entire FastA record and return the sequence header/sequence'''
        header = ''
        sequence = ''
        
        with self.doOpen() as fileH:
            
            header = ''
            sequence = ''
            
            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>') :
                line = fileH.readline()
            header = line[1:].rstrip()

            for line in fileH:
                if line.startswith ('>'):
                    yield header,sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else :
                    sequence += ''.join(line.rstrip().split()).upper()

        yield header,sequence


class NucParams:
    """
    This class provides methods to calculate various properties of a nucleotide sequence.
    """
    rnaCodonTable = {
    # RNA codon table
    # U
    'UUU': 'F', 'UCU': 'S', 'UAU': 'Y', 'UGU': 'C',  # UxU
    'UUC': 'F', 'UCC': 'S', 'UAC': 'Y', 'UGC': 'C',  # UxC
    'UUA': 'L', 'UCA': 'S', 'UAA': '-', 'UGA': '-',  # UxA
    'UUG': 'L', 'UCG': 'S', 'UAG': '-', 'UGG': 'W',  # UxG
    # C
    'CUU': 'L', 'CCU': 'P', 'CAU': 'H', 'CGU': 'R',  # CxU
    'CUC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',  # CxC
    'CUA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',  # CxA
    'CUG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',  # CxG
    # A
    'AUU': 'I', 'ACU': 'T', 'AAU': 'N', 'AGU': 'S',  # AxU
    'AUC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',  # AxC
    'AUA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',  # AxA
    'AUG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',  # AxG
    # G
    'GUU': 'V', 'GCU': 'A', 'GAU': 'D', 'GGU': 'G',  # GxU
    'GUC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',  # GxC
    'GUA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',  # GxA
    'GUG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'  # GxG
    }
    dnaCodonTable = {key.replace('U','T'):value for key, value in rnaCodonTable.items()}

    def __init__ (self, inString=''):
        """
        Initializes the NucParams object with a nucleotide sequence and calculates the composition of nucleotides,
        codons, and amino acids.
        """
        self.codonComp = {}
        for codon in NucParams.rnaCodonTable.keys():
            self.codonComp[codon] = 0
        self.aaComp = {}
        for aa in NucParams.rnaCodonTable.values():
            self.aaComp[aa] = 0
        self.nucComp = {}
        for nuc in ["A", "C", "G", "T", "U", "N"]:
            self.nucComp[nuc] = 0
        
        self.addSequence(inString)
    
    def addSequence (self, inSeq):
        """
        Adds a nucleotide sequence to the object and updates the composition of nucleotides, codons, and amino acids.
        """
        seq = inSeq.upper().replace("T", "U")   # Convert all input seq to RNA 
        for i in range(0,len(seq),3):
            tempCodon = seq[i:i+3]   #Generate codon strings of 3 nucleotides/letters
            if tempCodon in NucParams.rnaCodonTable.keys():   # Check if string is a valid codon present in rnaCodonTable
                self.codonComp[tempCodon] += 1
                aaVal = NucParams.rnaCodonTable[tempCodon]
                self.aaComp[aaVal] += 1    #increment count
        
        for letter in seq:   
            if letter in self.nucComp.keys():   # Check if nucleotide, increment count
                self.nucComp[letter] += 1
            
    def aaComposition(self):
        """
        Returns the composition of amino acids in the nucleotide sequence.
        """
        return self.aaComp
    
    def nucComposition(self):
        """
        Returns the composition of nucleotides in the nucleotide sequence.
        """
        return self.nucComp
    
    def codonComposition(self):
        """
        Returns the composition of codons in the nucleotide sequence.
        """
        return self.codonComp
    
    def nucCount(self):
        """
        Returns the total count of nucleotides in the nucleotide sequence.
        """
        return sum(self.nucComp.values())
    

""" x = NucParams()
x.addSequence("AUGAUGAUGAUGACG")
print(x.aaComposition())
print(x.codonComposition())
print(x.nucComposition()) """


import numpy
class ProteinParam :
    """
    A class that calculates various properties of a protein sequence including number of amino acids, molecular weight,
    molar Extinction coefficient, mass Extinction coefficient, theoretical P.I., and Amino acid composition."
    """
# These tables are for calculating:
#     molecular weight (aa2mw), along with the mol. weight of H2O (mwH2O)
#     absorbance at 280 nm (aa2abs280)
#     pKa of positively charged Amino Acids (aa2chargePos)
#     pKa of negatively charged Amino acids (aa2chargeNeg)
#     and the constants aaNterm and aaCterm for pKa of the respective termini
#  Feel free to move these to appropriate methods as you like

# As written, these are accessed as class attributes, for example:
# ProteinParam.aa2mw['A'] or ProteinParam.mwH2O

    aa2mw = {
        'A': 89.093,  'G': 75.067,  'M': 149.211, 'S': 105.093, 'C': 121.158,
        'H': 155.155, 'N': 132.118, 'T': 119.119, 'D': 133.103, 'I': 131.173,
        'P': 115.131, 'V': 117.146, 'E': 147.129, 'K': 146.188, 'Q': 146.145,
        'W': 204.225,  'F': 165.189, 'L': 131.173, 'R': 174.201, 'Y': 181.189
        }

    mwH2O = 18.015
    aa2abs280= {'Y':1490, 'W': 5500, 'C': 125}

    aa2chargePos = {'K': 10.5, 'R':12.4, 'H':6}
    aa2chargeNeg = {'D': 3.86, 'E': 4.25, 'C': 8.33, 'Y': 10}
    aaNterm = 9.69
    aaCterm = 2.34

    def __init__ (self, protein):
        """
        Initialize the ProteinParam class with a protein sequence.
        Sets the protein sequence and initializes the compositionAA dictionary.
        """
        self.protein_sequence = protein
        self.compositionAA = {'A': 0, 'C' : 0, 'D' : 0, 'E': 0, 'F': 0, 'G': 0, 'H': 0, 'I': 0, 'L': 0, 'K': 0, 
                         'M': 0, 'N': 0, 'P': 0, 'Q': 0, 'R': 0, 'S': 0, 'T': 0, 'V': 0, 'Y': 0, 'W': 0}
        for i in self.compositionAA:
            if i in self.protein_sequence:
                self.compositionAA[i] = self.protein_sequence.count(i)  #populate dictionary self.compositionAA with protein and count
    

    def aaCount (self):
        """
        Count the number of valid amino acids in the protein sequence.
        """
        counter = 0
        for letter in self.protein_sequence:
            if letter in self.compositionAA:
                counter += 1
        return counter

    def pI (self):
        """
        Calculate and return the isoelectric point (pI) of the protein sequence.
        """
        iso = numpy.inf
        for ph in range(0, 1401):
            phVal = ph / 100
            if abs(self._charge_(phVal)) < iso:    #search for ph closest to 0
                iso = self._charge_(phVal)
                bestpH = phVal
        return bestpH
    
    def aaComposition (self) :
        """
        Calculate and return the composition of amino acids in the protein sequence.
        """
        return self.compositionAA

    def _charge_ (self, pH):
        """
        Calculate and return the net charge of the protein sequence at a given pH.
        """
        poschargeSum = 0
        negchargeSum = 0
        
        for aa, pKa in ProteinParam.aa2chargePos.items() : 
            poschargeSum += (self.compositionAA[aa] * 10**pKa) / (10**pKa + 10**pH)
        
        for aa, pKa in ProteinParam.aa2chargeNeg.items() : 
            negchargeSum += (self.compositionAA[aa] * 10**pH) / (10**pKa + 10**pH)
       
        nTermcharge = 10**ProteinParam.aaNterm / (10**ProteinParam.aaNterm + 10**pH) #Assume there is only one N-terminus
        cTermcharge = 10**pH / (10**ProteinParam.aaCterm + 10**pH) # Assume there is only one C-terminus
        poschargeSum += nTermcharge
        negchargeSum += cTermcharge
        return poschargeSum  - negchargeSum
        
    def molarExtinction (self):
        """
        Calculate and return the molar extinction coefficient of the protein sequence.
        """
        molarE = 0
        for aa in ProteinParam.aa2abs280:
            molarE += ProteinParam.aa2abs280[aa] * self.compositionAA[aa]
        return molarE

    def massExtinction (self):
        """
        Calculate and return the mass extinction coefficient of the protein sequence.
        """
        myMW =  self.molecularWeight()
        return self.molarExtinction() / myMW if myMW else 0.0

    def molecularWeight (self):
        """
        Calculate and return the molecular weight of the protein sequence.
        """
        aamolecularweightSum = 0
        for aa in self.protein_sequence:
            if aa in self.compositionAA:
                aamolecularweightSum += ProteinParam.aa2mw[aa]
        watermwSum = self.aaCount() * ProteinParam.mwH2O
        molecWeight = (aamolecularweightSum - watermwSum) + ProteinParam.mwH2O
        return molecWeight

# Please do not modify any of the following.  This will produce a standard output that can be parsed
    
