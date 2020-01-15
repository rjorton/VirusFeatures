package cpb_machine;

import java.io.*;
import java.util.HashMap;

public class CPB_Machine {

    public final static int AMINO_NUM=21;
    
    public static String[][] aaCode, aminoAcids, aminoCodons;
    public static int[] aminoCodonCounts, codCounts, aaCounts;
    public static int[][] codPairCounts, aaPairCounts;
    public static HashMap<String,Integer> aaLookUp, codLookUp;

    public static String seq="", seqName="", prevName="";
    public static int seqCount=0, lineCount=0;
    public static int totLength=0, totSeqs=0, totCP=0;
    
    public static int a=0, c=0, g=0, t=0, n=0;
    public static int totCod=0, nullCod=0, stops=0;
    public static double gc=0, at=0;
    
    public static int aaN=0, acN=0, agN=0, atN=0;
    public static int caN=0, ccN=0, cgN=0, ctN=0;
    public static int gaN=0, gcN=0, ggN=0, gtN=0;
    public static int taN=0, tcN=0, tgN=0, ttN=0;

    public static double apa=0, apc=0, apg=0, apt=0;
    public static double cpa=0, cpc=0, cpg=0, cpt=0;
    public static double gpa=0, gpc=0, gpg=0, gpt=0;
    public static double tpa=0, tpc=0, tpg=0, tpt=0;
    
    public static int aBr=0, cBr=0, gBr=0, tBr=0, nBr=0;
    public static int aNonBr=0, cNonBr=0, gNonBr=0, tNonBr=0, nNonBr=0;
    
    public static int aaBrN=0, acBrN=0, agBrN=0, atBrN=0;
    public static int caBrN=0, ccBrN=0, cgBrN=0, ctBrN=0;
    public static int gaBrN=0, gcBrN=0, ggBrN=0, gtBrN=0;
    public static int taBrN=0, tcBrN=0, tgBrN=0, ttBrN=0;
    public static int aaNonBrN=0, acNonBrN=0, agNonBrN=0, atNonBrN=0;
    public static int caNonBrN=0, ccNonBrN=0, cgNonBrN=0, ctNonBrN=0;
    public static int gaNonBrN=0, gcNonBrN=0, ggNonBrN=0, gtNonBrN=0;
    public static int taNonBrN=0, tcNonBrN=0, tgNonBrN=0, ttNonBrN=0;
    
    public static double apaBr=0, apcBr=0, apgBr=0, aptBr=0;
    public static double cpaBr=0, cpcBr=0, cpgBr=0, cptBr=0;
    public static double gpaBr=0, gpcBr=0, gpgBr=0, gptBr=0;
    public static double tpaBr=0, tpcBr=0, tpgBr=0, tptBr=0;
    public static double apaNonBr=0, apcNonBr=0, apgNonBr=0, aptNonBr=0;
    public static double cpaNonBr=0, cpcNonBr=0, cpgNonBr=0, cptNonBr=0;
    public static double gpaNonBr=0, gpcNonBr=0, gpgNonBr=0, gptNonBr=0;
    public static double tpaNonBr=0, tpcNonBr=0, tpgNonBr=0, tptNonBr=0;
    
    public static int totBr=0, totNonBr=0, totBrN=0, totNonBrN=0;
   
    public static double codBias[], aaBias[], cpb[][];
    public static double cpbAv=0, cpbSum=0, cpbMin=0, trueCpbMin=0;
    public static int cpbCount=0;
    
    public static String filename="", outFilename="";
    
    public static String selFilename="";
    public static String selData[][];
    
    //these aren't used in analysis anymore (but still outputted, so left to preserve column format)
    public static int bad=0, notGood=0;
  
    
    public static void main(String[] args) {

        //arg1 is the FASTA sequences - single line format ESSENTIAL
        //arg2 is the metadata file - accessions being the link - updated so can run without it
        
        if(args.length==1) {
            filename=args[0];
        }
        else if(args.length==2) {
            filename=args[0];
            selFilename=args[1];
        }
        else {
            System.out.println("Error - incorrect parameters "+args.length);
            System.out.println("java -jar CPB_Machine.jar seqs.fasta");
            System.out.println("java -jar CPB_Machine.jar seqs.fasta virusFile.txt");
            System.exit(0);
        }

        if(filename.contains("."))
            outFilename=filename.substring(0, filename.lastIndexOf("."))+"_dat.txt";
        else
            outFilename=filename+"_dat.txt";
        
        System.out.println("In Seq File= "+filename);
        System.out.println("Out Dat File = "+outFilename);
        
        //Read in the number of entries in the selFile - list of accessions with associated metadata to link to
        //This data is not needed for any of the seq analysis - it was requested to be added to link the data to metadata to ease some downstream stuff
        if(args.length==2) {
            File inFile = new File(selFilename);
            lineCount=0;
            try {
                BufferedReader input =  new BufferedReader(new FileReader(inFile));

                try {
                    String line = null;

                    while (( line = input.readLine()) != null) {
                        lineCount++;
                    }
                  }
                finally {
                    input.close();
                }
            }
            catch (IOException ex) {
                ex.printStackTrace();
            }
            selData=new String[lineCount][4];
            System.out.println(lineCount+" entries in "+selFilename);
            lineCount=0;
            try {
                BufferedReader input =  new BufferedReader(new FileReader(inFile));

                try {
                    String line = null;

                    while (( line = input.readLine()) != null) {
                        String splits[]=line.split("\t");
                        selData[lineCount][0]=splits[0];//accession
                        selData[lineCount][1]=splits[1];//taxID
                        selData[lineCount][2]=splits[2];//virus name
                        selData[lineCount][3]=splits[3];//family
                        lineCount++;
                    }
                  }
                finally {
                    input.close();
                }
            }
            catch (IOException ex) {
                ex.printStackTrace();
            }
        }
        else {
            //create a dummy one entry array if file not provided
            selData=new String[1][4];
            selData[0][0]="BLAHortBLAH";
            selData[0][1]="BLAHortBLAH";
            selData[0][2]="BLAHortBLAH";
            selData[0][3]="BLAHortBLAH";
        }
        lineCount=0;
        
        createCode();
        
        codCounts=new int[64];//the counts of each codon
        codBias=new double[64];//codon bias - for an AA with multiple codons - how biased is each one (essentially simple prop)
        
        aaCounts=new int[21];//the counts of each AA
        aaBias=new double[21];//for all the AA in a ORF - how biased is each one (essentially simple prop)
        
        codPairCounts=new int[64][64];//the counts of codon-codon pairs
        cpb=new double[64][64];//codon pair bias
        aaPairCounts=new int[21][21];//the counts of AA-AA pairs (needed for cpb calculations)
        
        //Open the output file - as outputting data for each seq as its processed
        //Originally this was built for 100Ks of sequences so didn't want to read in all to memory
        try {
            FileWriter fstream = new FileWriter(outFilename);
            BufferedWriter out = new BufferedWriter(fstream);
        
            //write the output file header
            out.write(genHeader());
            
            try {
                File seqFile = new File(filename);
                BufferedReader input =  new BufferedReader(new FileReader(seqFile));

                try {
                    String line = null;

                    seqName="BLAHSTART";
                    lineCount=0;

                    while (( line = input.readLine()) != null) {
                        
                        if(line.indexOf(">")==0) {
                            
                            //added this in so will run on any seq without throwing error 
                            //as opposed to seqs only with >lcl| and _cds_
                            int from=0;
                            int to=line.length();
                            
                            //Seqs are in form >lcl|GenomeAccession_cds_blahblah
                            //We need to strip out the GenomeAccession to link to metadata
                            if(line.indexOf(">lcl|")==0 & line.indexOf("_cds_")>0) {
                                from=5;
                                to=line.indexOf("_cds_");
                            }
                            /*
                            //used for separate processing
                            //>Genome1_Seq1
                            //>Genome1_Seq2
                            //>Genome2_Seq1
                            //>Genome2_Seq2
                            //etc
                            else if(line.indexOf("_")>0) {
                                from=0;
                                to=line.indexOf("_");
                            }
                            */
                            
                            //if the sequence assession is the same as the prev seq - they will be merged together as same genome
                            prevName=seqName;
                            seqName=line.substring(from,to);
                            
                            if(lineCount==0)
                                prevName=seqName;
                            
                            seqCount++;
                        }
                        else {   
                            //if the current sequence name/accesion does not equal the previous
                            //means new sequence - output the existing the data and then clear it, before re-analysing next seq
                            if(!seqName.equals(prevName)) {
                                out.write(genOutput());
                                clearDat();
                            }
                            
                            //previously kicked out (and reported) any sequence with any non ACGT characters
                            //any ambiguities will now be counted as N in analysis
                            //added the below replaces in for distri - the seqs I used were checked for ambi characters before
                            seq=line.toUpperCase().replace("U", "T").replace("-", "");
                                
                            //used to be quite strict with amiguities etc
                            //left multiple of 3 in as good way to flag up errorenous seqs (incorrect slippage etc)
                            if(seq.length()%3!=0)
                                System.out.println("Error - seq length not divisible by 3 "+seqName+" "+seq.length()+" "+line);
                            else
                                analyseSeq(seq);
                        }
                        
                        lineCount++;
                    }
                }
                finally {
                    input.close();
                }

            }
            catch (IOException ex) {
                ex.printStackTrace();
            }
        
            //catch the last sequence in the file - which will be missed in above loop
            prevName=seqName;
            out.write(genOutput());
            out.close();
        }

        catch (Exception e) {
            e.printStackTrace();
            System.err.println("Error: " + e.getMessage());
        }
        
    }//end of main()
    
    public static void clearDat() {
        
        //nucleotide counts
        a=c=g=t=n=0;
        //gc and at content [not dinucs]
        gc=at=0;
        //dinuc counts
        aaN=acN=agN=atN=caN=ccN=cgN=ctN=gaN=gcN=ggN=gtN=taN=tcN=tgN=ttN=0;
        //dinuc freqs
        apa=apc=apg=apt=cpa=cpc=cpg=cpt=gpa=gpc=gpg=gpt=tpa=tpc=tpg=tpt=0;
        //bridge nucleotide counts
        aBr=cBr=gBr=tBr=nBr=0;
        //non-bridge nucleotide counts
        aNonBr=cNonBr=gNonBr=tNonBr=nNonBr=0;
        //bridge dinuc counts
        aaBrN=acBrN=agBrN=atBrN=0;
        caBrN=ccBrN=cgBrN=ctBrN=0;
        gaBrN=gcBrN=ggBrN=gtBrN=0;
        taBrN=tcBrN=tgBrN=ttBrN=0;
        //non-bridge dinuc counts
        aaNonBrN=acNonBrN=agNonBrN=atNonBrN=0;
        caNonBrN=ccNonBrN=cgNonBrN=ctNonBrN=0;
        gaNonBrN=gcNonBrN=ggNonBrN=gtNonBrN=0;
        taNonBrN=tcNonBrN=tgNonBrN=ttNonBrN=0;
        //bridge dinuc freqs
        apaBr=apcBr=apgBr=aptBr=0;
        cpaBr=cpcBr=cpgBr=cptBr=0;
        gpaBr=gpcBr=gpgBr=gptBr=0;
        tpaBr=tpcBr=tpgBr=tptBr=0;
        //non-bridge dinuc freqs
        apaNonBr=apcNonBr=apgNonBr=aptNonBr=0;
        cpaNonBr=cpcNonBr=cpgNonBr=cptNonBr=0;
        gpaNonBr=gpcNonBr=gpgNonBr=gptNonBr=0;
        tpaNonBr=tpcNonBr=tpgNonBr=tptNonBr=0;
        //total bridge & non-bridge dinucs, and nucleotides(N)
        totBr=totNonBr=totBrN=totNonBrN=0;
        //total codons, null codons (with non ACGT), total codon-pairs, read through stops
        totCod=nullCod=totCP=stops=0;
        //cpb sum, count, average
        cpbSum=cpbCount=0;
        cpbAv=0;
        //total sequence length (sum of all seqs on genome) and total number of seqs in genome
        totLength=totSeqs=0;
        //cpbMin and truCpbMin [as -9999 used for missing - see later]
        cpbMin=trueCpbMin=1000000;
        
        //codon counts
        for(int i=0;i<codCounts.length;i++)
            codCounts[i]=0;
        
        //codon biases - for each codon - simple freq bias as to how much it is used to code for its AA
        for(int i=0;i<codBias.length;i++)
            codBias[i]=0;
        
        //aa counts
        for(int i=0;i<aaCounts.length;i++)
            aaCounts[i]=0;
        
        //aa biases - simple freq out of AA count over total AA in genome
        for(int i=0;i<aaBias.length;i++)
            aaBias[i]=0;
        
        //aa piar counts - needed for cpb
        for(int i=0;i<aaPairCounts.length;i++)
            for(int j=0;j<aaPairCounts.length;j++)
                aaPairCounts[i][j]=0;
        
        //codon pair counts - needed for cpb
        for(int i=0;i<codPairCounts.length;i++)
            for(int j=0;j<codPairCounts[i].length;j++)
                codPairCounts[i][j]=0;
        
        //codon pair bias scores
        for(int i=0;i<cpb.length;i++)
            for(int j=0;j<cpb.length;j++)
                cpb[i][j]=0;

    }//end of clearDat()
    
    
    public static String genHeader() {
        
        String head="TaxID\tSpecies\tSeqName\t";
        head+="Good\tComplete\t";
        head+="Seqs\tSeqLength\tCodons\tBadCodons\tCodonPairs\tStops\t";
        head+="A\tC\tG\tT\tN\t";
        head+="GC\tAT\tCpG\tUpA\t";

        for(int j=0;j<aaBias.length;j++)
            head+=aminoAcids[0][j]+"-Bias\t";

        for(int j=0;j<codBias.length;j++)
            head+=aaCode[0][j]+"-Bias\t";

        for(int j=0;j<codPairCounts.length;j++)
            for(int k=0;k<codPairCounts[j].length;k++)
                    head+=aaCode[0][j]+"("+aaCode[1][j]+")-"+aaCode[0][k]+"("+aaCode[1][k]+")\t";
   
        //technically yhe CpsMin could be /2 if >0 - not used anymore but left in [after -9999 introduced]
        head+="CpsMin*2\tCpsAv";
        //standard dinucs - without cpg and upa as they are reported near start of columns
        //UpT to UpU change
        head+="\tApA\tApC\tApG\tApU\tCpA\tCpC\tCpU\tGpA\tGpC\tGpG\tGpU\tUpC\tUpG\tUpU";
        //bridge dinucs
        head+="\tbrApA\tbrApC\tbrApG\tbrApU\tbrCpA\tbrCpC\tbrCpG\tbrCpU\tbrGpA\tbrGpC\tbrGpG\tbrGpU\tbrUpA\tbrUpC\tbrUpG\tbrUpU";
        //non-bridge dinucs
        head+="\tNonBrApA\tNonBrApC\tNonBrApG\tNonBrApU\tNonBrCpA\tNonBrCpC\tNonBrCpG\tNonBrCpU\tNonBrGpA\tNonBrGpC\tNonBrGpG\tNonBrGpU\tNonBrUpA\tNonBrUpC\tNonBrUpG\tNonBrUpU";
        
        head+=System.lineSeparator();
        
        return head;
        
    }//end of genHeader()
    
    
    public static String genOutput() {
        
        //dat is the string with all the data to be returned for output
        String dat="";
        
        //first check if this sequence is in the metadata (selData) via accession number - if yes output the taxID[1] and virus name[2]
        boolean metaCheck=false;
        
        for(int i=0;i<selData.length;i++) {
            if(selData[i][0].equals(prevName)) {
               dat+=selData[i][1]+"\t"+selData[i][2]+"\t"; 
               metaCheck=true;
               break;
            }
        }
        
        //if not in metadata just output question marks
        if(!metaCheck) {
            dat+="?\t?\t";
        }
        
        //bad and notGood are relics from the past (polymerases etc) - left in to preserve column format if column numbering is essential downstream
        dat+=prevName+"\t"+bad+"\t"+notGood+"\t";
        dat+=totSeqs+"\t"+totLength+"\t"+totCod+"\t"+nullCod+"\t"+totCP+"\t"+" "+stops+"\t";
        
        //nucleotide freqs
        dat+=((double)a/(double)totLength)+"\t";
        dat+=((double)c/(double)totLength)+"\t";
        dat+=((double)g/(double)totLength)+"\t";
        dat+=((double)t/(double)totLength)+"\t";
        dat+=((double)n/(double)totLength)+"\t";
        
        //gc/at content
        dat+=gc+"\t"+at+"\t";
        //cpg & upa[tpa] - only did these two originally
        dat+=cpg+"\t"+tpa+"\t";

        //AA biases
        for(int j=0;j<aaBias.length;j++)
            dat+=aaBias[j]+"\t";

        //Codon biases
        for(int j=0;j<codBias.length;j++)
            dat+=codBias[j]+"\t";
    
        //CPB scores
        for(int j=0;j<cpb.length;j++) {
            for(int k=0;k<cpb[j].length;k++) {
                dat+=cpb[j][k]+"\t";
            }
        }
        
        //CPB Min and av
        dat+=cpbMin+"\t"+cpbAv+"\t";
        
        //dinucelotide obs/exp [except cpg and tpa]
        dat+=apa+"\t"+apc+"\t"+apg+"\t"+apt+"\t";
        dat+=cpa+"\t"+cpc+"\t"+cpt+"\t";//removed cpg - as included above
        dat+=gpa+"\t"+gpc+"\t"+gpg+"\t"+gpt+"\t";
        dat+=tpc+"\t"+tpg+"\t"+tpt+"\t";//removed tpa - as included above
        
        //bridge dinucs
        dat+=apaBr+"\t"+apcBr+"\t"+apgBr+"\t"+aptBr+"\t";
        dat+=cpaBr+"\t"+cpcBr+"\t"+cpgBr+"\t"+cptBr+"\t";
        dat+=gpaBr+"\t"+gpcBr+"\t"+gpgBr+"\t"+gptBr+"\t";
        dat+=tpaBr+"\t"+tpcBr+"\t"+tpgBr+"\t"+tptBr+"\t";
        
        //NonBridge dinucs
        dat+=apaNonBr+"\t"+apcNonBr+"\t"+apgNonBr+"\t"+aptNonBr+"\t";
        dat+=cpaNonBr+"\t"+cpcNonBr+"\t"+cpgNonBr+"\t"+cptNonBr+"\t";
        dat+=gpaNonBr+"\t"+gpcNonBr+"\t"+gpgNonBr+"\t"+gptNonBr+"\t";
        dat+=tpaNonBr+"\t"+tpcNonBr+"\t"+tpgNonBr+"\t"+tptNonBr;//removed last tab
        
        dat+=System.lineSeparator();//changed from /n for platform indy
                
        return dat;

    }//end of genOutput()
    
    
    public static void analyseSeq(String seq) {

        //as coding seqs of same genome are merged together keep a running total of all the data
        //and calculate the biases as don't know if this will be the last sequence
        //originally being used for 100Ks of sequences so didn't want to read evetyhing into memory
        
        totLength+=seq.length();
        totSeqs++;
        
        //Nucleotides
        for(int i=0;i<seq.length();i++) {
            
            //nucleotide counts
            if(seq.charAt(i)=='A')
                a++;
            else if(seq.charAt(i)=='C')
                c++;
            else if(seq.charAt(i)=='G')
                g++;
            else if(seq.charAt(i)=='T')
                t++;
            else
                n++;
            
            //bridge nucleotide counts - bridge = 2-3,5-6,8-9 etc
            if(i%3==0 & i>0) { 
                
                totBr++;
                if(seq.charAt(i-1)=='A')
                    aBr++;
                else if(seq.charAt(i-1)=='C')
                    cBr++;
                else if(seq.charAt(i-1)=='G')
                    gBr++;
                else if(seq.charAt(i-1)=='T')
                    tBr++;
                else
                    nBr++;
                
                totBr++;
                if(seq.charAt(i)=='A')
                    aBr++;
                else if(seq.charAt(i)=='C')
                    cBr++;
                else if(seq.charAt(i)=='G')
                    gBr++;
                else if(seq.charAt(i)=='T')
                    tBr++;
                else
                    nBr++;
            }
        }
        
        //Non-bridge actually includes all nucls 123-123, bridge is 3-1, non bridge is 1-2,2-3,1-2,2-3
        aNonBr=a;
        cNonBr=c;
        gNonBr=g;
        tNonBr=t;
        nNonBr=n;
        
        //totLength should be same as sum[a,c,g,t,n] so totNonBr = totLength
        totNonBr=totLength;
        
        //standard gc/at content - nothing to do with dinucs
        gc=(double)(g+c)/(double)totLength;
        at=(double)(a+t)/(double)totLength;
        
        //Count the dinucleotides - start at 1 as doing current and previous base
        for(int i=1;i<seq.length();i++) {
            if(seq.charAt(i-1)=='A' & seq.charAt(i)=='A')
                aaN++;
            else if(seq.charAt(i-1)=='A' & seq.charAt(i)=='C')
                acN++;
            else if(seq.charAt(i-1)=='A' & seq.charAt(i)=='G')
                agN++;
            else if(seq.charAt(i-1)=='A' & seq.charAt(i)=='T')
                atN++;
            else if(seq.charAt(i-1)=='C' & seq.charAt(i)=='A')
                caN++;
            else if(seq.charAt(i-1)=='C' & seq.charAt(i)=='C')
                ccN++;
            else if(seq.charAt(i-1)=='C' & seq.charAt(i)=='G')
                cgN++;
            else if(seq.charAt(i-1)=='C' & seq.charAt(i)=='T')
                ctN++;
            else if(seq.charAt(i-1)=='G' & seq.charAt(i)=='A')
                gaN++;
            else if(seq.charAt(i-1)=='G' & seq.charAt(i)=='C')
                gcN++;
            else if(seq.charAt(i-1)=='G' & seq.charAt(i)=='G')
                ggN++;
            else if(seq.charAt(i-1)=='G' & seq.charAt(i)=='T')
                gtN++;
            else if(seq.charAt(i-1)=='T' & seq.charAt(i)=='A')
                taN++;
            else if(seq.charAt(i-1)=='T' & seq.charAt(i)=='C')
                tcN++;
            else if(seq.charAt(i-1)=='T' & seq.charAt(i)=='G')
                tgN++;
            else if(seq.charAt(i-1)=='T' & seq.charAt(i)=='T')
                ttN++;
            
            //BRIDGE dinculs, so 3,6,9 etc are codon pos 1, we pull out 3-1 dinucl
            if(i%3==0) {
                totBrN++;
                if(seq.charAt(i-1)=='A' & seq.charAt(i)=='A')
                    aaBrN++;
                else if(seq.charAt(i-1)=='A' & seq.charAt(i)=='C')
                    acBrN++;
                else if(seq.charAt(i-1)=='A' & seq.charAt(i)=='G')
                    agBrN++;
                else if(seq.charAt(i-1)=='A' & seq.charAt(i)=='T')
                    atBrN++;
                else if(seq.charAt(i-1)=='C' & seq.charAt(i)=='A')
                    caBrN++;
                else if(seq.charAt(i-1)=='C' & seq.charAt(i)=='C')
                    ccBrN++;
                else if(seq.charAt(i-1)=='C' & seq.charAt(i)=='G')
                    cgBrN++;
                else if(seq.charAt(i-1)=='C' & seq.charAt(i)=='T')
                    ctBrN++;
                else if(seq.charAt(i-1)=='G' & seq.charAt(i)=='A')
                    gaBrN++;
                else if(seq.charAt(i-1)=='G' & seq.charAt(i)=='C')
                    gcBrN++;
                else if(seq.charAt(i-1)=='G' & seq.charAt(i)=='G')
                    ggBrN++;
                else if(seq.charAt(i-1)=='G' & seq.charAt(i)=='T')
                    gtBrN++;
                else if(seq.charAt(i-1)=='T' & seq.charAt(i)=='A')
                    taBrN++;
                else if(seq.charAt(i-1)=='T' & seq.charAt(i)=='C')
                    tcBrN++;
                else if(seq.charAt(i-1)=='T' & seq.charAt(i)=='G')
                    tgBrN++;
                else if(seq.charAt(i-1)=='T' & seq.charAt(i)=='T')
                    ttBrN++;
            }
        }
        
        //Nonbridge (1-2,2-3) is total (1-2,2-3,3-1) minus bridge (3-1)
        aaNonBrN=aaN-aaBrN;
        acNonBrN=acN-acBrN;
        agNonBrN=agN-agBrN;
        atNonBrN=atN-atBrN;
        
        caNonBrN=caN-caBrN;
        ccNonBrN=ccN-ccBrN;
        cgNonBrN=cgN-cgBrN;
        ctNonBrN=ctN-ctBrN;
        
        gaNonBrN=gaN-gaBrN;
        gcNonBrN=gcN-gcBrN;
        ggNonBrN=ggN-ggBrN;
        gtNonBrN=gtN-gtBrN;
        
        taNonBrN=taN-taBrN;
        tcNonBrN=tcN-tcBrN;
        tgNonBrN=tgN-tgBrN;
        ttNonBrN=ttN-ttBrN;
        
        //totLength-totSeqs is the total dinucls (length-1) for each seq
        //tot nonbridge dinculs is total dinucs minus bridge dinucs
        totNonBrN=totLength-totSeqs-totBrN;
        
        //calculate dinucl biases
        //these are the observed/expected
        //add option tp output/evaluate raw dinucs as a measure
        apa=((double)aaN/((double)totLength-(double)totSeqs))/(((double)a/(double)totLength)*((double)a/(double)totLength));
        apc=((double)acN/((double)totLength-(double)totSeqs))/(((double)a/(double)totLength)*((double)c/(double)totLength));
        apg=((double)agN/((double)totLength-(double)totSeqs))/(((double)a/(double)totLength)*((double)g/(double)totLength));
        apt=((double)atN/((double)totLength-(double)totSeqs))/(((double)a/(double)totLength)*((double)t/(double)totLength));
        
        cpa=((double)caN/((double)totLength-(double)totSeqs))/(((double)c/(double)totLength)*((double)a/(double)totLength));
        cpc=((double)ccN/((double)totLength-(double)totSeqs))/(((double)c/(double)totLength)*((double)c/(double)totLength));
        cpg=((double)cgN/((double)totLength-(double)totSeqs))/(((double)c/(double)totLength)*((double)g/(double)totLength));
        cpt=((double)ctN/((double)totLength-(double)totSeqs))/(((double)c/(double)totLength)*((double)t/(double)totLength));
        
        gpa=((double)gaN/((double)totLength-(double)totSeqs))/(((double)g/(double)totLength)*((double)a/(double)totLength));
        gpc=((double)gcN/((double)totLength-(double)totSeqs))/(((double)g/(double)totLength)*((double)c/(double)totLength));
        gpg=((double)ggN/((double)totLength-(double)totSeqs))/(((double)g/(double)totLength)*((double)g/(double)totLength));
        gpt=((double)gtN/((double)totLength-(double)totSeqs))/(((double)g/(double)totLength)*((double)t/(double)totLength));
        
        tpa=((double)taN/((double)totLength-(double)totSeqs))/(((double)t/(double)totLength)*((double)a/(double)totLength));
        tpc=((double)tcN/((double)totLength-(double)totSeqs))/(((double)t/(double)totLength)*((double)c/(double)totLength));
        tpg=((double)tgN/((double)totLength-(double)totSeqs))/(((double)t/(double)totLength)*((double)g/(double)totLength));
        tpt=((double)ttN/((double)totLength-(double)totSeqs))/(((double)t/(double)totLength)*((double)t/(double)totLength));
        
        //bridge
        //observed expected of the bridge, given the nucleotides at the bridge
        //again, idea to evalue plain old raw freqs of dinucs over total dinucs
        apaBr=((double)aaBrN/((double)totBrN))/(((double)aBr/(double)totBr)*((double)aBr/(double)totBr));
        apcBr=((double)acBrN/((double)totBrN))/(((double)aBr/(double)totBr)*((double)cBr/(double)totBr));
        apgBr=((double)agBrN/((double)totBrN))/(((double)aBr/(double)totBr)*((double)gBr/(double)totBr));
        aptBr=((double)atBrN/((double)totBrN))/(((double)aBr/(double)totBr)*((double)tBr/(double)totBr));
        
        cpaBr=((double)caBrN/((double)totBrN))/(((double)cBr/(double)totBr)*((double)aBr/(double)totBr));
        cpcBr=((double)ccBrN/((double)totBrN))/(((double)cBr/(double)totBr)*((double)cBr/(double)totBr));
        cpgBr=((double)cgBrN/((double)totBrN))/(((double)cBr/(double)totBr)*((double)gBr/(double)totBr));
        cptBr=((double)ctBrN/((double)totBrN))/(((double)cBr/(double)totBr)*((double)tBr/(double)totBr));
        
        gpaBr=((double)gaBrN/((double)totBrN))/(((double)gBr/(double)totBr)*((double)aBr/(double)totBr));
        gpcBr=((double)gcBrN/((double)totBrN))/(((double)gBr/(double)totBr)*((double)cBr/(double)totBr));
        gpgBr=((double)ggBrN/((double)totBrN))/(((double)gBr/(double)totBr)*((double)gBr/(double)totBr));
        gptBr=((double)gtBrN/((double)totBrN))/(((double)gBr/(double)totBr)*((double)tBr/(double)totBr));
        
        tpaBr=((double)taBrN/((double)totBrN))/(((double)tBr/(double)totBr)*((double)aBr/(double)totBr));
        tpcBr=((double)tcBrN/((double)totBrN))/(((double)tBr/(double)totBr)*((double)cBr/(double)totBr));
        tpgBr=((double)tgBrN/((double)totBrN))/(((double)tBr/(double)totBr)*((double)gBr/(double)totBr));
        tptBr=((double)ttBrN/((double)totBrN))/(((double)tBr/(double)totBr)*((double)tBr/(double)totBr));
        
        //NonBridge
        apaNonBr=((double)aaNonBrN/((double)totNonBrN))/(((double)aNonBr/(double)totNonBr)*((double)aNonBr/(double)totNonBr));
        apcNonBr=((double)acNonBrN/((double)totNonBrN))/(((double)aNonBr/(double)totNonBr)*((double)cNonBr/(double)totNonBr));
        apgNonBr=((double)agNonBrN/((double)totNonBrN))/(((double)aNonBr/(double)totNonBr)*((double)gNonBr/(double)totNonBr));
        aptNonBr=((double)atNonBrN/((double)totNonBrN))/(((double)aNonBr/(double)totNonBr)*((double)tNonBr/(double)totNonBr));
        
        cpaNonBr=((double)caNonBrN/((double)totNonBrN))/(((double)cNonBr/(double)totNonBr)*((double)aNonBr/(double)totNonBr));
        cpcNonBr=((double)ccNonBrN/((double)totNonBrN))/(((double)cNonBr/(double)totNonBr)*((double)cNonBr/(double)totNonBr));
        cpgNonBr=((double)cgNonBrN/((double)totNonBrN))/(((double)cNonBr/(double)totNonBr)*((double)gNonBr/(double)totNonBr));
        cptNonBr=((double)ctNonBrN/((double)totNonBrN))/(((double)cNonBr/(double)totNonBr)*((double)tNonBr/(double)totNonBr));
        
        gpaNonBr=((double)gaNonBrN/((double)totNonBrN))/(((double)gNonBr/(double)totNonBr)*((double)aNonBr/(double)totNonBr));
        gpcNonBr=((double)gcNonBrN/((double)totNonBrN))/(((double)gNonBr/(double)totNonBr)*((double)cNonBr/(double)totNonBr));
        gpgNonBr=((double)ggNonBrN/((double)totNonBrN))/(((double)gNonBr/(double)totNonBr)*((double)gNonBr/(double)totNonBr));
        gptNonBr=((double)gtNonBrN/((double)totNonBrN))/(((double)gNonBr/(double)totNonBr)*((double)tNonBr/(double)totNonBr));
        
        tpaNonBr=((double)taNonBrN/((double)totNonBrN))/(((double)tNonBr/(double)totNonBr)*((double)aNonBr/(double)totNonBr));
        tpcNonBr=((double)tcNonBrN/((double)totNonBrN))/(((double)tNonBr/(double)totNonBr)*((double)cNonBr/(double)totNonBr));
        tpgNonBr=((double)tgNonBrN/((double)totNonBrN))/(((double)tNonBr/(double)totNonBr)*((double)gNonBr/(double)totNonBr));
        tptNonBr=((double)ttNonBrN/((double)totNonBrN))/(((double)tNonBr/(double)totNonBr)*((double)tNonBr/(double)totNonBr));
                
        //Count codons and AAs
        for(int i=0;i<seq.length();i+=3) {
            //in case a non-multiple of 3 sequence comes through, stop before last
            if((i+3)>seq.length())
                break;
            
            String cod=""+seq.charAt(i)+seq.charAt(i+1)+seq.charAt(i+2);
            
            if(codLookUp.get(cod)==null)
                nullCod++;
            else {
                codCounts[codLookUp.get(cod)]++;
                aaCounts[aaLookUp.get(aaCode[1][codLookUp.get(cod)])]++;
            }
            
            totCod++;
        }
        
        //AA bias - for each AA - simplge freq over all AA (total codons)
        for(int i=0;i<aaCounts.length;i++) {
            aaBias[i]=(double)aaCounts[i]/(double)totCod;
        }
        
        //Codon Bias - for each codon, frequency for how it used to code for its AA (considering all codons of that AA)
        for(int i=0;i<codCounts.length;i++) {
            if(codCounts[i]>0 & aaCounts[aaLookUp.get(aaCode[1][i])]==0)
                System.out.println("Error "+i+" codCounts>0 but AA=0 codCount="+codCounts[i]+" "+aaLookUp.get(aaCode[1][i]));
            
            if(aaCounts[aaLookUp.get(aaCode[1][i])]>0)
                codBias[i]=(double)codCounts[i]/(double)aaCounts[aaLookUp.get(aaCode[1][i])];
        }
        
        //Codon Pair & AA Pair counts
        for(int i=3;i<seq.length();i+=3) {
            if((i+3)>seq.length())
                break;
            
            String cod1=""+seq.charAt(i-3)+seq.charAt(i-2)+seq.charAt(i-1);
            String cod2=""+seq.charAt(i)+seq.charAt(i+1)+seq.charAt(i+2);
            
            //keep count of readthrough stops
            if(cod1.equals("TAG")|cod1.equals("TAA")|cod1.equals("TGA"))
                stops++;
            
            //ignore codon pairs with Ns or ambiguitites
            if(codLookUp.get(cod1)!=null & codLookUp.get(cod2)!=null) {
                codPairCounts[codLookUp.get(cod1)][codLookUp.get(cod2)]++;
                aaPairCounts[aaLookUp.get(aaCode[1][codLookUp.get(cod1)])][aaLookUp.get(aaCode[1][codLookUp.get(cod2)])]++;
            }
            
            totCP++;
        }
                
        //Codon Pair Bias coleman formula
        //CPB = ln [codonPairCount / [[(codon1Count x codon2Count) / (aa1Count * aa2Count)] x aaPairCount]]
        for(int i=0;i<cpb.length;i++) {
            for(int j=0;j<cpb[i].length;j++) {
                if(codPairCounts[i][j]>0) {
                    cpb[i][j]=Math.log((double)codPairCounts[i][j]/((((double)codCounts[i]*(double)codCounts[j])/(aaCounts[aaLookUp.get(aaCode[1][i])]*aaCounts[aaLookUp.get(aaCode[1][j])]))*aaPairCounts[aaLookUp.get(aaCode[1][i])][aaLookUp.get(aaCode[1][j])]));
                    //Keep a running total of all CPBs across all pairs and across all coding seqs
                    cpbSum+=cpb[i][j];
                    cpbCount++;
                }
            }
        }
        
        //V1 empty codonPairCounts/CPBs were simply left alone so 0
        //V2 changed to cpbMin*2 [before -9999 introdcued]
        //V3 changed to cpbAv for codonPairs where the AApair is not observed (primarily those involving stops)
        //and -9999 for when codonPair is 0 but AApair>0
        //Thinking - should cases where AApair = 0 be penalised even more, rather than some average, would -9999 or -9999*2 actually be better - this would also act as a anchor/link between viruses - the average is arbitary and won't link between viruses
        //In terms of codon pair bias, they are not biased [as they are not used]
        //But in terms of AA pair bias they are (if the AA are observed in the seq that is)
        //Maybe AApair bias could be a measure
        //The change from V1 to V2 to V3 was not based on any results, thought to a better way rather than having a 0 (which would imply no bias)
        cpbAv=cpbSum/(double)cpbCount;
        
        //we don't use the cpbMin/trueCpbMin anymore [after the whole -9999] - but left in
        for(int i=0;i<cpb.length;i++) {
            for(int j=0;j<cpb[i].length;j++) {
                if(cpb[i][j]<trueCpbMin & cpb[i][j]!=cpbMin)
                    trueCpbMin=cpb[i][j];
            }
        }
        
        if(trueCpbMin<0)
            cpbMin=trueCpbMin*2;
        else
            cpbMin=trueCpbMin/2;
                
        for(int i=0;i<cpb.length;i++) {
            for(int j=0;j<cpb[i].length;j++) {
                if(codPairCounts[i][j]==0) {
                    //if codon pair count is 0, and the AA pair count is also zero - set to the cpbAv
                    //mainly StopCodon-NonStopCodon - read through stops could be removed as rare
                    if(aaPairCounts[aaLookUp.get(aaCode[1][i])][aaLookUp.get(aaCode[1][j])]==0) {
                        cpb[i][j]=cpbAv;
                    }
                    //if codon pair count is 0, but there are AA pair counts - set to -9999
                    else {
                        cpb[i][j]=-9999;
                    }
                }
            }
        }
        
    }//end of analyseSeq
    
    
    public static void createCode() {
         
        aminoAcids=new String[2][AMINO_NUM];//stores the AA codes and names
        aminoAcids[0][0]="L";
        aminoAcids[1][0]="leu";
        
        aminoAcids[0][1]="P";
        aminoAcids[1][1]="pro";
        
        aminoAcids[0][2]="H";
        aminoAcids[1][2]="his";
        
        aminoAcids[0][3]="Q";
        aminoAcids[1][3]="gln";
        
        aminoAcids[0][4]="R";
        aminoAcids[1][4]="arg";
        
        aminoAcids[0][5]="I";
        aminoAcids[1][5]="ile";

        aminoAcids[0][6]="M";
        aminoAcids[1][6]="met";//start
        
        aminoAcids[0][7]="T";
        aminoAcids[1][7]="thr";
        
        aminoAcids[0][8]="N";
        aminoAcids[1][8]="asn";
        
        aminoAcids[0][9]="K";
        aminoAcids[1][9]="lys";
        
        aminoAcids[0][10]="S";
        aminoAcids[1][10]="ser";
      
        aminoAcids[0][11]="V";
        aminoAcids[1][11]="val";
        
        aminoAcids[0][12]="A";
        aminoAcids[1][12]="ala";
        
        aminoAcids[0][13]="D";
        aminoAcids[1][13]="asp";
        
        aminoAcids[0][14]="E";
        aminoAcids[1][14]="glu";
        
        aminoAcids[0][15]="G";
        aminoAcids[1][15]="gly";
        
        aminoAcids[0][16]="F";
        aminoAcids[1][16]="phe";
        
        aminoAcids[0][17]="Y";
        aminoAcids[1][17]="tyr";

        aminoAcids[0][18]="C";
        aminoAcids[1][18]="cys";
        
        aminoAcids[0][19]="W";
        aminoAcids[1][19]="trp";
        
        aminoAcids[0][20]="X";
        aminoAcids[1][20]="stp";//stop
        
        //check for duplicate AAs - as list above was manually defined
        for(int i=0;i<aminoAcids[0].length;i++) {
            for(int j=i+1;j<aminoAcids[0].length;j++) {
                if(aminoAcids[0][i].equalsIgnoreCase(aminoAcids[0][j])) {
                    System.out.println("Duplicate AAs - "+i+ " "+j+" - "+aminoAcids[0][i]+" and "+aminoAcids[0][j]);
                }
                
                if(aminoAcids[1][i].equalsIgnoreCase(aminoAcids[1][j])) {
                    System.out.println("Duplicate AAs - "+i+ " "+j+" - "+aminoAcids[1][i]+" and "+aminoAcids[1][j]);
                }
            }
        }
        
        aaCode=new String[2][64];
        
        aaCode[0][0]="AAA";
        aaCode[1][0]="K";
        aaCode[0][1]="AAC";
        aaCode[1][1]="N";
        aaCode[0][2]="AAG";
        aaCode[1][2]="K";
        aaCode[0][3]="AAT";
        aaCode[1][3]="N";
        
        aaCode[0][4]="ACA";
        aaCode[1][4]="T";
        aaCode[0][5]="ACC";
        aaCode[1][5]="T";
        aaCode[0][6]="ACG";
        aaCode[1][6]="T";
        aaCode[0][7]="ACT";
        aaCode[1][7]="T";
        
        aaCode[0][8]="AGA";
        aaCode[1][8]="R";
        aaCode[0][9]="AGC";
        aaCode[1][9]="S";
        aaCode[0][10]="AGG";
        aaCode[1][10]="R";
        aaCode[0][11]="AGT";
        aaCode[1][11]="S";
        
        aaCode[0][12]="ATA";
        aaCode[1][12]="I";
        aaCode[0][13]="ATC";
        aaCode[1][13]="I";
        aaCode[0][14]="ATG";
        aaCode[1][14]="M";//start
        aaCode[0][15]="ATT";
        aaCode[1][15]="I";
        
        aaCode[0][16]="CAA";
        aaCode[1][16]="Q";
        aaCode[0][17]="CAC";
        aaCode[1][17]="H";
        aaCode[0][18]="CAG";
        aaCode[1][18]="Q";
        aaCode[0][19]="CAT";
        aaCode[1][19]="H";
        
        aaCode[0][20]="CCA";
        aaCode[1][20]="P";
        aaCode[0][21]="CCC";
        aaCode[1][21]="P";
        aaCode[0][22]="CCG";
        aaCode[1][22]="P";
        aaCode[0][23]="CCT";
        aaCode[1][23]="P";
        
        aaCode[0][24]="CGA";
        aaCode[1][24]="R";
        aaCode[0][25]="CGC";
        aaCode[1][25]="R";
        aaCode[0][26]="CGG";
        aaCode[1][26]="R";
        aaCode[0][27]="CGT";
        aaCode[1][27]="R";
        
        aaCode[0][28]="CTA";
        aaCode[1][28]="L";
        aaCode[0][29]="CTC";
        aaCode[1][29]="L";
        aaCode[0][30]="CTG";
        aaCode[1][30]="L";
        aaCode[0][31]="CTT";
        aaCode[1][31]="L";
        
        aaCode[0][32]="GAA";
        aaCode[1][32]="E";
        aaCode[0][33]="GAC";
        aaCode[1][33]="D";
        aaCode[0][34]="GAG";
        aaCode[1][34]="E";
        aaCode[0][35]="GAT";
        aaCode[1][35]="D";
        
        aaCode[0][36]="GCA";
        aaCode[1][36]="A";
        aaCode[0][37]="GCC";
        aaCode[1][37]="A";
        aaCode[0][38]="GCG";
        aaCode[1][38]="A";
        aaCode[0][39]="GCT";
        aaCode[1][39]="A";
        
        aaCode[0][40]="GGA";
        aaCode[1][40]="G";
        aaCode[0][41]="GGC";
        aaCode[1][41]="G";
        aaCode[0][42]="GGG";
        aaCode[1][42]="G";
        aaCode[0][43]="GGT";
        aaCode[1][43]="G";
        
        aaCode[0][44]="GTA";
        aaCode[1][44]="V";
        aaCode[0][45]="GTC";
        aaCode[1][45]="V";
        aaCode[0][46]="GTG";
        aaCode[1][46]="V";
        aaCode[0][47]="GTT";
        aaCode[1][47]="V";
        
        aaCode[0][48]="TAA";
        aaCode[1][48]="X";//stop
        aaCode[0][49]="TAC";
        aaCode[1][49]="Y";
        aaCode[0][50]="TAG";
        aaCode[1][50]="X";//stop
        aaCode[0][51]="TAT";
        aaCode[1][51]="Y";
        
        aaCode[0][52]="TCA";
        aaCode[1][52]="S";
        aaCode[0][53]="TCC";
        aaCode[1][53]="S";
        aaCode[0][54]="TCG";
        aaCode[1][54]="S";
        aaCode[0][55]="TCT";
        aaCode[1][55]="S";
        
        aaCode[0][56]="TGA";
        aaCode[1][56]="X";//stop
        aaCode[0][57]="TGC";
        aaCode[1][57]="C";
        aaCode[0][58]="TGG";
        aaCode[1][58]="W";//tryp - only codon for tryp
        aaCode[0][59]="TGT";
        aaCode[1][59]="C";
        
        aaCode[0][60]="TTA";
        aaCode[1][60]="L";
        aaCode[0][61]="TTC";
        aaCode[1][61]="F";
        aaCode[0][62]="TTG";
        aaCode[1][62]="L";
        aaCode[0][63]="TTT";
        aaCode[1][63]="F";
        
        //as above was manually created - some checks for sanity
        for(int i=0;i<aaCode[0].length;i++) {
            
            //check for blanks
            if((aaCode[0][i].equalsIgnoreCase(""))) {
                System.out.println("Blank codon - "+i+" - "+aaCode[0][i]+" "+aaCode[1][i]);
            }
            
            boolean aaTest=false;
            //check that the AA symbols used in aaCode are present in aminoAcids[][]
            for(int j=0;j<aminoAcids[0].length;j++) {
                if(aaCode[1][i].equalsIgnoreCase(aminoAcids[0][j])) {
                    aaTest=true;
                    break;
                }
            }
            if(!aaTest)
                System.out.println("Codon AA not found in aminoAcids- "+i+" - "+aaCode[0][i]+" "+aaCode[1][i]);
            
            //check codons are unique
            for(int j=i+1;j<aaCode[0].length;j++) {
                if(aaCode[0][i].equalsIgnoreCase(aaCode[0][j])) {
                    System.out.println("Duplicate codons - "+i+" "+j+" - "+aaCode[0][i]+" "+aaCode[0][j]);
                }
            }
        }
        
        //check each aminoAcid[][] appears in codons aaCode[][]
        for(int i=0;i<aminoAcids[0].length;i++) {
            boolean aaTest=false;
            
            for(int j=0;j<aaCode[1].length;j++) {
                if(aminoAcids[0][i].equalsIgnoreCase(aaCode[1][j])) {
                    aaTest=true;
                    break;
                }
            }
            
            if(!aaTest) {
                System.out.println("aminoAcid not found in codons - "+i+" - "+aminoAcids[0][i]+" "+aminoAcids[1][i]);
            }
        }
        
        //create AA index lookup - where are they located in the arrays
        aaLookUp=new HashMap<String,Integer>();
        for(int i=0;i<aminoAcids[0].length;i++) {
            aaLookUp.put(aminoAcids[0][i], i);
        }
        
        //create codon index lookup - where are they located in the arrays
        codLookUp=new HashMap<String,Integer>();
        for(int i=0;i<aaCode[0].length;i++) {
            codLookUp.put(aaCode[0][i], i);
        }
        
        //For each AA - see how many codons code it - and populate their data
        aminoCodons=new String[AMINO_NUM][10];//ten is arbitrary and excess - maximum number of codons that an individual AA has
        //blank the data out first
        for(int i=0;i<aminoCodons.length;i++) {
            for(int j=0;j<aminoCodons[i].length;j++) {
                aminoCodons[i][j]="";
            }
        }
        aminoCodonCounts=new int[AMINO_NUM];//stores the number of codons each AA has
        for(int i=0;i<aminoAcids[0].length;i++) {

            int thisCount=0;
            
            for(int j=0;j<aaCode[1].length;j++) {
                if(aminoAcids[0][i].equalsIgnoreCase(aaCode[1][j])) {
                    aminoCodons[i][thisCount]=aaCode[0][j];
                    thisCount++;
                }
            }
            
            aminoCodonCounts[i]=thisCount;
        }
 
    }//end of CreateCode()   
    
}
