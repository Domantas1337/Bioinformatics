import os
from Bio import SeqIO, AlignIO
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Align.Applications import ClustalwCommandline
from Bio.Align import substitution_matrices
from Bio.Align.Applications import MafftCommandline


def readFastaFile(filepath):
    with open(filepath, 'r') as file:
        sequences = [record for record in SeqIO.parse(file, 'fasta')]
    return sequences


def getHumanSequenceFromMsa(alignment, human_id = "NP_000468.1"):
    for record in alignment:
        if record.id == human_id:
            return record
    return None 

def getAlignedSequences(fileName):
    return AlignIO.read(fileName, "clustal")

def calculateSimilarity(fragment, secondFragment):
    blosum80 = substitution_matrices.load("BLOSUM80")
    score = 0

    for firstAcid, secondAcid in zip(fragment, secondFragment):
        if (firstAcid, secondAcid) in blosum80:
            score += blosum80[(firstAcid, secondAcid)]
        elif (secondAcid, firstAcid) in blosum80:
            score += blosum80[(secondAcid, firstAcid)]
    return score



def findUniqueSequence(humanSequence, alignment):
    uniqueFragment = None
    similarFragment = None
    lowestSimilarity = float('inf')
    largestSimilarity = float('-inf')

    for i in range(len(humanSequence.seq) - 14):
        fragment = humanSequence.seq[i:i+15]
        totalSimilarity = 0

        for record in alignment:
            if record.id != humanSequence.id:
                totalSimilarity += calculateSimilarity(fragment, record.seq[i:i+15])

        if totalSimilarity < lowestSimilarity:
            uniqueFragment = fragment
            lowestSimilarity = totalSimilarity
        if totalSimilarity > largestSimilarity:
            similarFragment = fragment
            largestSimilarity = totalSimilarity

    return uniqueFragment, lowestSimilarity, similarFragment, largestSimilarity

def alignSequences(sequences, outputFile):
    with open(outputFile, "w") as f:
        for sequence in sequences:
            f.write(f">{sequence.id}\n{str(sequence.seq)}\n")

    useMafft = MafftCommandline("mafft.bat", input = outputFile)
    stdout, _ = useMafft()

    outputFileForAlignment = "alignedSequences.fasta"
    with open(outputFileForAlignment, "w") as aligned:
        aligned.write(stdout)
    
    try:
        alignment = AlignIO.read(outputFileForAlignment, "fasta")
    except ValueError as ve:
        return None

    return alignment


def main():
    humanFastaFile = 'albuminas.fasta'
    blastFastaFile = 'blast_results.fasta'
    alignedFastaFile = 'aligned_sequences.fasta'

    

    humanProteinSequence = readFastaFile(humanFastaFile)
    comparableProteinSequences = readFastaFile(blastFastaFile)
    
    allSequences = humanProteinSequence + comparableProteinSequences

    alignSequences(allSequences, alignedFastaFile)

    aligned_sequences = list(SeqIO.parse("alignedSequences.fasta", "fasta"))

    result = findUniqueSequence(aligned_sequences[0], aligned_sequences)
    print(result)



if __name__ == "__main__":
    main()