Program to create various features from viral coding sequences, ranging from nucleotide, dinucleotide and AA frequencies, to dinucleotide, codon and codon pair biases.

This is a java program that is run as:

java -jar CPB_Machine.jar InputSequenceFile.fasta VirusAccessions.txt

The CPB_Machine program will create a file called InputSequenceFile_dat.txt with all the features in it.

The input sequence file can have numerous sequences in it, but the program expects sequences to be in FASTA single line format (i.e. one header line followed by one sequence line), sequences must be in frame coding sequences as the program is calculating codon and AA frequencies.
The VirusAccesions.txt has a list of virus accesion numbers, taxID, and Species names (this is currently being deleted).

The program currently looks for Genbank style headers in order to link coding sequences from the same viral genome/segment together, taking the above example it takes the substring from character number 5 to the start of the _cds_ to identify the genbank accession number. Sequences from the same accession are analysed (and outputted) together:

lcl|NC_001437.1_cds_NP_059434.1_1 description

lcl|NC_001437.1_cds_YP_006355435.1_2 description

An example input sequence file (picrona_new.fasta) and example program output file (picorna_new_dat.txt) are provided.
