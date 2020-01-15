### VirusFeatures

This is in the process of being replaced by [SeqFeats](https://github.com/rjorton/SeqFeats) in python.

Program to create various features from viral coding sequences, ranging from nucleotide, dinucleotide and AA frequencies, to dinucleotide, codon and codon pair biases.

This is a java program that is run as:

---
    java -jar CPB_Machine.jar InputSequenceFile.fasta
    or
    java -jar CPB_Machine.jar InputSequenceFile.fasta VirusAccessions.txt
---

The VirusAccesions.txt file was just used to add a virus name and taxID to the output, it is not needed to run. The CPB_Machine program will create a file called InputSequenceFile_dat.txt with all the features in it.

The input sequence file can have numerous sequences in it, but the program expects sequences to be in FASTA single line format (i.e. one header line followed by one sequence line), there should be no blank lines, sequences must be in frame coding sequences as the program is calculating codon and AA frequencies. It will throw an error if sequences are not a multiple of 3.

The VirusAccesions.txt has a list of virus accesion numbers, taxID, and Species names (this is currently being deleted), if you accession is not present in the file, a '?' symbol should be outputted in the virus name column.

The program currently looks for Genbank style headers in order to link coding sequences from the same viral genome/segment together, taking the above example it takes the substring from character number 5 (after the pipe | symbol) to the start of the _cds_ to identify the genbank accession number. Sequences from the same accession are analysed (and outputted) together:

lcl|NC_001437.1_cds_NP_059434.1_1 description

lcl|NC_001437.1_cds_YP_006355435.1_2 description

If an lcl| and _cds_ is not present in the header, it will use the whole sequence header. An example input sequence file (picrona_new.fasta) and example program output file (picorna_new_dat.txt) are provided.

Fields outputted:

* TaxID, Species (obtained from the VirusFeatures.txt file via accesion lookup)
* SeqName
* Good, Complete - these can be ignored
* Seqs - the number of sequences assigned to this virus
* SeqLength - total length of all sequences assigned to this virus
* Codons - total number of codons
* BadCodons - any codon with a N or amiguity code
* CodonPairs - total number of codon pairs
* Stops - total number of read through stops
* A, C, G, T, N - the freq of A, C, G, T and N in the sequences
* GC, AT - the GC and AT content of the sequences
* CpG, UpA - the observed/expected CpG and UpA ratios
* L-Bias, P-Bias ... X-Bias - the frequency of each AA in the sequence - X represents stops
* AAA-Bias, AAC-Bias ... TTT-Bias - the codon bias for each codon, simple frequency: number of times codon is used, divided by number of times corresponding AA is used
* AAA(K)-AAA(K) ... TTT(F)-TTT(F) - the codon pair bias using the formula of Coleman et al (2008) (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2754401/) 
* CpsMin2 - the minimum codon pair CPS score x 2
* CpsAv - the average CPS across all codon pairs
* ApA ... UpU - the observed/expected dinucleotide biases
* brApA ... brUpU - the observed/expected dinucleotide biases at the bridge between codons
* NonBrApA ... NonBrUpU - the observed/expected dinucleotide biases at the non-bridge between codons



