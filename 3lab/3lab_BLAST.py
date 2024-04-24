import os
from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML

def readFastaFile(filepath):
    with open(filepath, 'r') as file:
        return SeqIO.read(file, 'fasta')

def filterBlastResults(blastResults, sequenceLength):
    
    blastOutput = NCBIXML.parse(blastResults)
    sequences = []

    for blastRecord in blastOutput:
        for alignment in blastRecord.alignments:
            for hsp in alignment.hsps:
                coverage = hsp.align_length / sequenceLength
                if coverage >= 0.8:
                    title = alignment.title.split()[0] 
                    sequences.append((title, hsp.sbjct))

    return sequences

def performRemoteBlast(sequence):
    results = NCBIWWW.qblast('blastp', 'swissprot', sequence.seq, 
                                   entrez_query='mammals[Organism]')
    return results

def saveSequences(sequences, output_file):
    with open(output_file, 'w') as file:
        for seq in sequences:
            file.write(f">{seq[0]}\n{seq[1]}\n")

def main():
    fastaFile = 'albuminas.fasta'
    outputFile = 'blast_results.fasta'

    # Read the FASTA file
    albuminSequence = readFastaFile(fastaFile)

    blastResults = performRemoteBlast(albuminSequence)

    sequenceLength = len(albuminSequence.seq)
    extractedSequences = filterBlastResults(blastResults, sequenceLength)

    saveSequences(extractedSequences, outputFile)


if __name__ == "__main__":
    main()