#!/usr/bin/env python3
# Name: Riya Tiloda (rtiloda)
# Group Members: Hema Prasanna, Akash Pandit

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
    
import sys
def main():
    """
    Main function to interact with the user and perform the calculations on the input protein sequence.
    """
    inString = input('protein sequence?')
    while inString :
        myParamMaker = ProteinParam(inString)
        myAAnumber = myParamMaker.aaCount()
        print ("Number of Amino Acids: {aaNum}".format(aaNum = myAAnumber))
        print ("Molecular Weight: {:.1f}".format(myParamMaker.molecularWeight()))
        print ("molar Extinction coefficient: {:.2f}".format(myParamMaker.molarExtinction()))
        print ("mass Extinction coefficient: {:.2f}".format(myParamMaker.massExtinction()))
        print ("Theoretical pI: {:.2f}".format(myParamMaker.pI()))
        #break
        print ("Amino acid composition:")
        
        if myAAnumber == 0 : myAAnumber = 1  # handles the case where no AA are present 
        
        for aa,n in sorted(myParamMaker.aaComposition().items(), 
                           key= lambda item:item[0]):
            print ("\t{} = {:.2%}".format(aa, n/myAAnumber))   #prints aa composition with correct formatting
    
        #inString = input('protein sequence?')
        break

if __name__ == "__main__":
    main()