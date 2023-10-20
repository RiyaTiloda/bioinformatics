from sequenceAnalysis import FastAreader as FastaReader
#Pseudo Code for OrfFinder:
#for loop that iterates through 3 frames (0 to 2):
#for loop that iterates through codons
#if start codon is encountered, add it to empty lists of start codons and dangling starts
#if a stop codon is encountered:
    #if there are no start codons in startCodonPositions and there are no ORFs in self.allORFs, it is a dangling stop. Add the orf to self.allORFs
    #if there are no start codons in startCodonPositions, but there are ORFs in self.allORFs, then it is a arbitrary stop codon
#if there is a start codon found before the stop codon, this segment of the sequence can be stored as an ORF. A list of the start codon pos., stop pos., and length is added to self.allORFs
#danglingStarts and startCodonPositions are set to empty lists after the ORF is added to self.allORFs
#if danglingStarts contains a start codon outside of the inner for loop, a dangling start has been found. The dangling start ORF is added to self.allORFs
#Now, we iterated through the reverse strand. We use the same nested for loop, but change the way in which the start and stop codon positions are added to self.allORFs.
#Since we are iterating through the reverse complement, start and stop codon positions must be stored according to the forward strand
#start position of ORF = (len(reverse sequence) - stop codon pos) + 1.
#end position of ORF = (len(reverse sequence) - start codon pos) + 1.

class OrfFinder:

    '''This class searches the seqeunce for ORFs. It identifies ORFs longer than or equal to 100 nucleotides.
    It identifies dangling stop and start ORFs. It contains a method which generates the reverse complement
    and searches for ORFs within the reverse complement, while storing positions in accordance to the first
    strand. ORF fram, start pos., end pos., length are stored in a nested list which is returned by the findORFs
    method.'''

    def __init__(self, seq):
        self.definedStartCodon = "ATG"
        self.definedStopCodons = ["TAG", "TAA", "TGA"]
        self.sequence = seq.replace(" ","")
        self.allORFs = []
        

    def seqReverse(self):
        '''Create reverse complement of first strand sequence.'''
        complement = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
        reverseComplement = ""
        for nucleotide in self.sequence: 
            reverseComplement += complement.get(nucleotide)
        return reverseComplement[::-1]
        

    def findORFs(self):
        '''Search for ORFs, store and return nested list self.allORFs'''
        
        danglingStarts = []
        reverseStarts = []
        stopCodonPosition = 0
        for frame in range(3):   #iterates through each frame

            startIndex = frame
            startCodonPositions = []

            #print("checkpoint 1")
            for nucIndex in range(startIndex, len(self.sequence), 3):   #iterates through sequence within a frame, by 3 to enabling access to each codon
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

                        if not self.allORFs:   #identify dangling stops
                            #print("no orfs", self.allORFs)
                            stopCodonPosition = nucIndex + 2
                            lengthofORF = stopCodonPosition
                            self.allORFs += [[frame + 1, 1, stopCodonPosition + 1, lengthofORF + 1]]
                            continue

                        else:   #arbitrary stop codon encountered
                            #print("list is empty")
                            continue

                    stopCodonPosition = nucIndex + 2
                    
                    #print("debug", startCodonPositions)

                    lengthofORF = stopCodonPosition - startCodonPositions[0]
                    self.allORFs += [[frame + 1, startCodonPositions[0] + 1, stopCodonPosition + 1, lengthofORF + 1]]
                    
                    startCodonPositions = [] #ORF was found, so reset startCodonPositions and danglingStarts
                    danglingStarts = []
                    #print("dangling starts cleared")
                
                if startCodonPositions:
                    danglingStarts += startCodonPositions
                    #print("Dangling:", danglingStarts)
                #If any start codons are left at this point where frame is over, then it is a dangling start
                
        #print("list of any remaining start positions", startCodonPositions)

        
            #print("Dangling s",danglingStarts)
            if danglingStarts: #account for dangling starts
                stopCodonPosition = len(self.sequence) - 1
                lengthofORF = stopCodonPosition - (danglingStarts[0] + 1)
                self.allORFs += [[danglingStarts[0]%3 + 1, startCodonPositions[0] + 1, stopCodonPosition + 1, lengthofORF + 1]]
            
        #REVERSE STRAND
        revSeq = self.seqReverse()   #store reverse strand
        revORFs = []
        #print("revseq:", revSeq)
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

                    #print("stop codon found", codon)
                    if not reverseStarts:
                        #print("no start codons found")

                        if not revORFs:   #no ORFs found yet, account for dangling stop
                            #print("no orfs", revORFs)
                            stopCodonPosition = nucIndex + 2
                            lengthofORF = stopCodonPosition
                            revORFs += [[-(frame + 1),  (len(revSeq)-stopCodonPosition), len(revSeq), lengthofORF + 1]]
                            continue

                        else:   #arbitrary stop codon found
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
            if danglingStarts:   #account for dangling starts
                stopCodonPosition = len(self.sequence)
                lengthofORF = stopCodonPosition - (danglingStarts[0] + 1)
                revORFs += [[-(danglingStarts[0]%3 + 1), 1, (len(revSeq) - danglingStarts[0]), lengthofORF + 1]]
            
        self.allORFs += revORFs

        return self.allORFs 
        

########################################################################
# CommandLine
########################################################################
class CommandLine() :
    '''
    Handle the command line, usage and help requests.

    CommandLine uses argparse, now standard in 2.7 and beyond. 
    it implements a standard command line argument parser with various argument options,
    a standard usage and help.

    attributes:
    all arguments received from the commandline using .add_argument will be
    avalable within the .args attribute of object instantiated from CommandLine.
    For example, if myCommandLine is an object of the class, and requiredbool was
    set as an option using add_argument, then myCommandLine.args.requiredbool will
    name that option.
 
    '''
    
    def __init__(self, inOpts=None) :
        '''
        Implement a parser to interpret the command line argv string using argparse.
        '''
        
        import argparse
        self.parser = argparse.ArgumentParser(description = 'Program prolog - a brief description of what this thing does', 
                                             epilog = 'Program epilog - some other stuff you feel compelled to say', 
                                             add_help = True, #default is True 
                                             prefix_chars = '-', 
                                             usage = '%(prog)s [options] -option1[default] <input >output'
                                             )
        self.parser.add_argument('-lG', '--longestGene', action = 'store', nargs='?', const=True, default=False, help='longest Gene in an ORF')
        self.parser.add_argument('-mG', '--minGene', type=int, choices= (100,200,300,500,1000), default=100, action = 'store', help='minimum Gene length')
        self.parser.add_argument('-s', '--start', action = 'append', default = ['ATG'],nargs='?', 
                                 help='start Codon') #allows multiple list options
        self.parser.add_argument('-t', '--stop', action = 'append', default = ['TAG','TGA','TAA'],nargs='?', help='stop Codon') #allows multiple list options
        self.parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')  
        if inOpts is None :
            self.args = self.parser.parse_args()
        else :
            self.args = self.parser.parse_args(inOpts)

########################################################################
# Main
########################################################################
   
def main(inFile = None, options = None):
    '''
    Create an instance of class OrfFinder, generate a nested list containing lists of the frame, start position, 
    stop position, and frame. Sort the ORFs in descending order of length, and print in the output file given by
    CommandLine. 
    '''
    thisCommandLine = CommandLine(options) #instantiate commandLine object with parameters(options)
    reader = FastaReader(inFile)
    
    for header, sequence in reader.readFasta():
        nucleotideSequence = OrfFinder(sequence)
        geneList = sorted(nucleotideSequence.findORFs(), key=lambda x:(-x[3], x[1])) #Sorted by descending length
        for orfList in geneList: 
            #print(orfList)
            if orfList[3] >= 100:
                print(f'{orfList[0]} {orfList[1]}..{orfList[2]}   {orfList[3]}')
    
if __name__ == "__main__":
    #delete this stuff if running from commandline
    #main(inFile = 'tass2.fa', options = ['-mG=300', '-lG'])    
    main()
