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
        seq = inSeq.upper().replace("T", "U")
        for i in range(0,len(seq),3):
            tempCodon = seq[i:i+3]
            if tempCodon in NucParams.rnaCodonTable.keys():
                self.codonComp[tempCodon] += 1
                aaVal = NucParams.rnaCodonTable[tempCodon]
                self.aaComp[aaVal] += 1
        
        for letter in seq:
            if letter in self.nucComp.keys():
                self.nucComp[letter] += 1
            
    def aaComposition(self):
        return self.aaComp
    
    def nucComposition(self):
        return self.nucComp
    
    def codonComposition(self):
        return self.codonComp
    
    def nucCount(self):
        return sum(self.nucComp.values())
    

""" x = NucParams()
x.addSequence("AUGAUGAUGAUGACG")
print(x.aaComposition())
print(x.codonComposition())
print(x.nucComposition()) """


#LAB05
#assume that sequence is taken as a string input
class OrfFinder:
    def __init__(self, seq):
        self.definedStartCodon = "ATG"
        self.definedStopCodons = ["TAG", "TAA", "TGA"]
        self.sequence = seq.replace(" ","")
        self.allORFs = []

    def seqReverse(self):
        """Create reverse complement of input DNA strand."""
        complement = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
        reverseComplement = ""
        for nucleotide in self.sequence: 
            reverseComplement += complement.get(nucleotide)
        return reverseComplement[::-1]
        

    def findORFs(self):
        
        #for loop to iterate through frames
        #frame = 0
        #startCodonPositions = []
        danglingStarts = []
        reverseStarts = []
        stopCodonPosition = 0
        for frame in range(3):

            startIndex = frame
            startCodonPositions = []

            #print("checkpoint 1")
            for nucIndex in range(startIndex, len(self.sequence), 3):
                #print("checkpoint 2")
                codon = self.sequence[nucIndex: nucIndex + 3]

                #print(codon)
                if codon == self.definedStartCodon:
                    #print("start codon found", codon)
                    startCodonPositions += [nucIndex]
                    #print("START CODON POSITIONS:", startCodonPositions)

                elif codon in self.definedStopCodons:

                    #print("stop codon found", codon)
                    if not startCodonPositions:

                        if not self.allORFs:   #no ORFs found yet
                            #print("no orfs", self.allORFs)
                            stopCodonPosition = nucIndex + 2
                            lengthofORF = stopCodonPosition
                            self.allORFs += [[frame + 1, 1, stopCodonPosition + 1, lengthofORF + 1]]
                            continue

                        else:
                            #print("list is empty")
                            continue

                    stopCodonPosition = nucIndex + 2
                    
                    #print("debug", startCodonPositions)

                    lengthofORF = stopCodonPosition - startCodonPositions[0]
                    self.allORFs += [[frame + 1, startCodonPositions[0] + 1, stopCodonPosition + 1, lengthofORF + 1]]
                    
                    startCodonPositions = []
                    danglingStarts = []
                    #print("dangling starts cleared")
                
                if startCodonPositions:
                    danglingStarts += startCodonPositions
                    #print("Dangling:", danglingStarts)
                #If any start codons are left at this point where frame is over, then it is a dangling start
                
        #print("list of any remaining start positions", startCodonPositions)

        
            #print("Dangling s",danglingStarts)
            if danglingStarts:
                stopCodonPosition = len(self.sequence) - 1
                lengthofORF = stopCodonPosition - (danglingStarts[0] + 1)
                self.allORFs += [[danglingStarts[0]%3 + 1, startCodonPositions[0] + 1, stopCodonPosition + 1, lengthofORF + 1]]
            
        #REVERSE STRAND
        revSeq = self.seqReverse()
        revORFs = []
        print("revseq:", revSeq)
        danglingStarts = []
        for frame in range(3):

            startIndex = frame
            reverseStarts = []

            #print("checkpoint 1")
            for nucIndex in range(startIndex, len(revSeq), 3):
                #print("checkpoint 2")
                codon = revSeq[nucIndex: nucIndex + 3]

                #print(codon)
                if codon == self.definedStartCodon:
                    #print("start codon found", codon)
                    reverseStarts += [nucIndex]
                    #print("START CODON POSITIONS:", reverseStarts)

                elif codon in self.definedStopCodons:

                    print("stop codon found", codon)
                    if not reverseStarts:
                        print("no start codons found")

                        if not revORFs:   #no ORFs found yet
                            print("no orfs", revORFs)
                            stopCodonPosition = nucIndex + 2
                            lengthofORF = stopCodonPosition
                            revORFs += [[-(frame + 1),  (len(revSeq)-stopCodonPosition), len(revSeq), lengthofORF + 1]]
                            continue

                        else:
                            #print("list is empty")
                            continue

                    stopCodonPosition = nucIndex + 2
                    
                    #print("debug", reverseStarts)

                    lengthofORF = stopCodonPosition - reverseStarts[0]
                    
                    start = len(revSeq) - stopCodonPosition
                    end = len(revSeq) - reverseStarts[0]
                    #print("START", start)
                    #print("STOP", end)

                    revORFs += [[-(frame + 1), start, end, lengthofORF + 1]]
                    
                    reverseStarts = []
                    danglingStarts = []
                
                if reverseStarts:
                    danglingStarts += reverseStarts
                #If any start codons are left at this point where frame is over, then it is a dangling start
                
        #print("list of any remaining start positions", reverseStarts)

        
            #print("Dangling s",danglingStarts)
            if danglingStarts:
                stopCodonPosition = len(self.sequence)
                lengthofORF = stopCodonPosition - (danglingStarts[0] + 1)
                revORFs += [[-(danglingStarts[0]%3 + 1), 1, (len(revSeq) - danglingStarts[0]), lengthofORF + 1]]
            
        self.allORFs += revORFs
        for sublist in self.allORFs:
            if sublist[3] < 100:
                self.allORFs.remove(sublist)

        return self.allORFs 
        
    