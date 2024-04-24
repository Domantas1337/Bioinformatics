import numpy as np
from scipy.signal import find_peaks
from Bio.Blast import NCBIWWW, NCBIXML
import matplotlib.pyplot as plt

RANGES = {
    'Sanger': (33, 73),
    'Solexa': (59, 104),
    'Illumina-1.3': (64, 104),
    'Illumina-1.5': (66, 105),
    'Illumina-1.8': (33, 74),
}

def getQualityScores(filePath):
    
    qualityScores = []
    index = 1

    with open(filePath, 'r') as file:
        
        for line in file:
            if index % 4 == 0:
                qualityScores.append(line.strip())

            index += 1

    return ''.join(qualityScores)

def findEncoding(qualityScores):

    minMaxSums = {}

    minValue = ord(min(qualityScores))
    maxValue = ord(max(qualityScores))
        
    for encoding, (asciiMin, asciiMax) in RANGES.items():
        minimumDifference = max(0, asciiMin - minValue)
        maximumDifference = max(0, maxValue - asciiMax)

        minMaxSums[encoding] = minimumDifference + maximumDifference
    
    return min(minMaxSums, key=minMaxSums.get)

def cgCount(sequence):
    
    numOfCG = 0

    for i in sequence:
        if i == 'C':
            numOfCG += 1
        if i == 'G':
            numOfCG += 1

    return (numOfCG / len(sequence)) * 100

def performRemoteBlast(sequence):
    
    result = NCBIWWW.qblast("blastn", "nt", sequence, hitlist_size=1, entrez_query="Bacteria [Organism]")
   
    blastRecord = NCBIXML.read(result)
    
    firstAlignment = blastRecord.alignments[0]
    title = firstAlignment.title

    organismName = title.split('|')[4].split(' ')[1:]
    organism = ' '.join(organismName)
    return organism

def findPeakIntervals(filePath):
    
    sequenceData = []

    with open(filePath, 'r') as file:
        while True:
            sequenceId = file.readline().strip()
            sequence = file.readline().strip()
            file.readline()  # Praleisti '+'
            file.readline()  # Praleisti quality score 

            if not sequenceId:
                break

            numberOfCGs = 0

            for base in sequence:
                if base == 'C':
                    numberOfCGs += 1
                if base == 'G':
                    numberOfCGs += 1

            sequenceData.append({"id": sequenceId[1:], "sequence": sequence, "percentage": (numberOfCGs / len(sequence)) * 100} )

    peakIntervals, peaks, binValues, frequencies = analyzeCGpercentages(sequenceData)

    return peaks, peakIntervals, binValues, frequencies, sequenceData


def findBlastSequences(peakIntervals, sequenceInformation):
    
    blastResults = []
    checkedSequences = []

    for peakInterval in peakIntervals:
    
        for information in sequenceInformation:
            sequencesInTheInterval = []
            
            for information in sequenceInformation:
                if information["percentage"] >= peakInterval[0] and information["percentage"] <= peakInterval[1]: 
                    sequencesInTheInterval.append((information["id"], information["sequence"]))


            for sequence in sequencesInTheInterval[:5]:
                if(not sequence[0] in checkedSequences):
                    organism = performRemoteBlast(sequence[1])
                    blastResults.append({'id': sequence[0],  'organism': organism})
                    checkedSequences.append(sequence[0])
                    print(str(sequence[0]) + '\n' + str(organism) + '\n')

    return blastResults

def analyzeCGpercentages(sequenceData):

    percentages = []
    
    for data in sequenceData:
        percentages.append(data["percentage"])

    frequencies , binValues = np.histogram(percentages, bins=25)
    peaks, _ = find_peaks(frequencies, prominence=5)
    
    peakIntervals = []

    for peak in peaks:
        peakIntervals.append((binValues[peak], binValues[peak + 1]))

    print("Peaks: " +  str(len(peaks)))

    return peakIntervals, peaks, binValues, frequencies

def plotGrpah(sequenceData, frequencies, binValues, peaks):

    percentages = []

    for data in sequenceData:
        percentages.append(data["percentage"])

    plt.hist(percentages, bins=25, alpha=0.7)

    ymin, ymax = plt.ylim()

    for peak in peaks:
        plt.plot(binValues[peak], frequencies[peak])

        offset = 10 + (frequencies[peak] / ymax) * 100
        plt.annotate(f"{binValues[peak]:.2f}%", 
                     (binValues[peak], frequencies[peak] + 20.5))
    
    plt.title('C/G nucleotide distribution', fontsize=18)
    plt.xlabel('Percentage of C/G nucleotides', fontsize=12)
    plt.ylabel('Read number', fontsize=12)

    plt.tight_layout()
    plt.savefig('my_plot.png')

def main():
    
    file = "reads_for_analysis.fastq"

    qualityScores = getQualityScores(file)
    encoding = findEncoding(qualityScores)
    print("Encoding: " + encoding)
    peaks, peakIntervals, binValues, frequencies, sequenceData = findPeakIntervals(file)
    print("Performing blast")
    
    #blastResults = findBlastSequences(peakIntervals, sequenceData)
    plotGrpah(sequenceData, frequencies, binValues, peaks)

    #print(blastResults)

if __name__ == "__main__":
    main()