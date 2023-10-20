def genDicts(targSet, setA, setB):
    dictTarget = {length : set() for length in targSet}
    dictA = {length : set() for length in setA}
    dictB = {length : set() for length in setB}

    for key_length in dictTarget:
        for orf in self.c1.ORFs:
            if orf[3] == key_length:
                seq = getORFseq(orf, ORFfinder.top)
                dictTarget[key_length].add(seq)

    for key_length in dictA:
        for orf in self.c1.ORFs:
            if orf[3] == key_length:
                seq = getORFseq(orf, ORFfinder.top)
                dictA[key_length].add(seq)

    for key_length in dictB:
        for orf in self.c1.ORFs:
            if orf[3] == key_length:
                seq = getORFseq(orf, ORFfinder.top)
                dictB[key_length].add(seq)

    return dictTarget, dictA, dictB



import matplotlib.pyplot as plt

# Sample data
data = [1, 2, 2, 3, 3, 3, 4, 4, 5, 5, 6]

# Create a histogram
plt.hist(data)

# Customize the histogram
plt.title("Histogram")
plt.xlabel("Values")
plt.ylabel("Frequency")

# Show the histogram
plt.show()

def genGraph(self, c1, c2):   #c1 c2 are the count of matches in fileA and fileB respectively to targFile ORFs
    dataSet1 = [c1, c2]
    plt.hist(dataSet1)
    plt.title("ORF Match frequency")
    plt.xlabel("File A Matches")
    plt.ylabel("File B Matches")

    


