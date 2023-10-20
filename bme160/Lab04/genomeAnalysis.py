#!/usr/bin/env python3
# Name: Riya Tiloda (rtiloda)
# Group Members: Hema Prasanna, Akash Pandit
#import ProteinParams
from sequenceAnalysis import FastAreader
from sequenceAnalysis import NucParams

def main (fileName=None):
    """
    Main function to analyze nucleotide sequences and calculate GC content and codon usage.
    
    Args:
        fileName (str): Name of the FASTA file to analyze (default is None).
    """
    
    myReader = FastAreader(fileName)   # Create a FastAreader object to read the FASTA file
    myNuc = NucParams()    # Create a NucParams object to store nucleotide and amino acid compositions
    for head, seq in myReader.readFasta() :   # Read each sequence from the FASTA file and add it to the NucParams object
        myNuc.addSequence(seq)
        
    # sort codons in alpha order, by Amino Acid
    print(f'sequence length = {myNuc.nucCount()/1000000:.2f} Mb')    # # Calculate and print the sequence length in Mb
    print()
    nucleotideDict = myNuc.nucComposition()
    gContent = nucleotideDict['G']
    cContent = nucleotideDict['C']
    #Calculate and print the GC content as a percentage
    if "G" in nucleotideDict and "C" in nucleotideDict:
        gcPercent = (gContent + cContent)/myNuc.nucCount()
    elif "G" in nucleotideDict:
        gcPercent = (gContent)/myNuc.nucCount()
    elif "C" in nucleotideDict:
        gcPercent = (cContent)/myNuc.nucCount()
    else:
        gcPercent = 0
        
    print(f'GC content = {gcPercent*100:.1f}%')
    print()
    
    
    # sort codons in alpha order, by Amino Acid
    codonCountDict = myNuc.codonComposition()
    aaCountDict = myNuc.aaComposition()
    sortedRNAtable = sorted(NucParams.rnaCodonTable.items(), key =lambda item: (item[1], item[0]))
    
    for codon, aa in sortedRNAtable:
        if codon in codonCountDict:
            print ('{:s} : {:s} {:5.1f} ({:6d})'.format(codon, aa, (codonCountDict[codon]/aaCountDict[aa])*100, codonCountDict[codon]))
    # calculate relative codon usage for each codon and print
    
if __name__ == "__main__":
    main('testGenome.fa') # make sure to change this in order to use stdin
    
    