/**
 * @(#)CallMpileup.java
 *
 * Copyright (c) 2009-2010 Daniel C. Koboldt and Washington University in St. Louis
 *
 * COPYRIGHT
 */

package net.sf.varscan;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.text.DecimalFormat;
import java.util.HashMap;
import java.lang.Math;

/**
 * A class for calling variants or consensus bases from a pileup file
 *
 * @version	2.3
 *
 * @author Daniel C. Koboldt <dkoboldt@genome.wustl.edu>
 *
 */
public class CallMpileup {

	public CallMpileup(String[] args, String callType)
	{
		// Define the usage message //
		String usage = "USAGE: java -jar VarScan.jar mpileup2cns [pileup file] OPTIONS\n" +
		"\tmpileup file - The SAMtools mpileup file\n" +
		"\n" +
		"\tOPTIONS:\n" +
		"\t--min-coverage\tMinimum read depth at a position to make a call [8]\n" +
		"\t--min-reads2\tMinimum supporting reads at a position to call variants [2]\n" +
		"\t--min-avg-qual\tMinimum base quality at a position to count a read [15]\n" +
		"\t--min-var-freq\tMinimum variant allele frequency threshold [0.01]\n" +
		"\t--min-freq-for-hom\tMinimum frequency to call homozygote [0.75]\n" +
		"\t--p-value\tDefault p-value threshold for calling variants [99e-02]\n" +
		"\t--strand-filter\tIgnore variants with >90% support on one strand [1]\n" +
		"\t--output-vcf\tIf set to 1, outputs in VCF format\n" +
		"\t--vcf-sample-list\tFor VCF output, a list of sample names in order, one per line\n" +
		"\t--variants\tReport only variant (SNP/indel) positions [0]";

		// Set parameter defaults //

		HashMap<String, String> params = VarScan.getParams(args);


		// Set up formatting for p-values //
		DecimalFormat pvalueFormat = new DecimalFormat("0.####E0");

		// If mpileup2snp or mpileup2indel was called, set the variants parameter //

		if(args[0].equals("mpileup2snp"))
			params.put("variants", "snp");

		if(args[0].equals("mpileup2indel"))
			params.put("variants", "indel");

		if(args[0].equals("mpileup2vcf"))
				params.put("output-vcf", "1");

		// Set parameter defaults //

		int minCoverage = 8;
		int minReads2 = 2;
		int minAvgQual = 15;
		double minVarFreq = 0.01;
		double minFreqForHom = 0.75;
		double pValueThreshold = 0.99;
		double strandPvalueThreshold = 0.01;
		boolean variantsOnly = false;
		boolean snpsOnly = false;
		boolean indelsOnly = false;
		boolean strandFilter = true;
		String sampleList = "";

		if(callType.equals("CNS"))
		{
			// Set more rigorous parameters for consensus calling
			minVarFreq = 0.20;
			pValueThreshold = 0.01;
		}

		// Adjust parameters based on user input //

		try
		{
			if(params.containsKey("min-coverage"))
				 minCoverage = Integer.parseInt(params.get("min-coverage"));

			if(params.containsKey("min-reads2"))
				 minReads2 = Integer.parseInt(params.get("min-reads2"));

			if(params.containsKey("min-var-freq"))
				 minVarFreq = Double.parseDouble(params.get("min-var-freq"));

			if(params.containsKey("min-freq-for-hom"))
				 minFreqForHom = Double.parseDouble(params.get("min-freq-for-hom"));

			if(params.containsKey("min-avg-qual"))
				 minAvgQual = Integer.parseInt(params.get("min-avg-qual"));

			if(params.containsKey("p-value"))
				 pValueThreshold = Double.parseDouble(params.get("p-value"));

			if(params.containsKey("strand-filter"))
			{
				int filter = Integer.parseInt(params.get("strand-filter"));
				if(filter > 0)
					strandFilter = true;
				else
					strandFilter = false;
			}

			if(params.containsKey("vcf-sample-list"))
			{
				File samplefile = new File(params.get("vcf-sample-list"));
				// Parse sample list //
				if(samplefile.exists())
	    		{
	    			BufferedReader in = new BufferedReader(new FileReader(samplefile));
	    			String line = "";
	    			if(in.ready())
	    			{
	    				while ((line = in.readLine()) != null)
	    				{
	    					String sampleName = line;
	    					if(sampleList.length() > 0)
	    						sampleList += "\t";
	    					sampleList += sampleName;
	    				}
	    			}
	    			else
	    			{
	    				System.err.println("Unable to open sample list");
	    			}

	    			in.close();
	    		}

				System.err.println("Got the following sample list: ");
				System.err.println(sampleList);
			}


			if(params.containsKey("variants"))
			{
				String variants = params.get("variants");

				// Determine type of variant reporting: all (default), SNPs, or indels //
				if(variants.equals("snp"))
				{
					snpsOnly = true;
					System.err.println("Only SNPs will be reported");
				}
				else if(variants.equals("indel"))
				{
					indelsOnly = true;
					System.err.println("Only indels will be reported");
				}
				else
				{
					variantsOnly = true;
					System.err.println("Only variants will be reported");
				}
			}

			if(params.containsKey("p-value"))
				 pValueThreshold = Double.parseDouble(params.get("p-value"));
			else
				System.err.println("Warning: No p-value threshold provided, so p-values will not be calculated");

			 System.err.println("Min coverage:\t" + minCoverage);
			 System.err.println("Min reads2:\t" + minReads2);
			 System.err.println("Min var freq:\t" + minVarFreq);
			 System.err.println("Min avg qual:\t" + minAvgQual);
			 System.err.println("P-value thresh:\t" + pValueThreshold);
		}
		catch(Exception e)
		{
	    	System.err.println("Input Parameter Threw Exception: " + e.getLocalizedMessage());
	    	e.printStackTrace(System.err);
	    	return;
		}

		// Print usage if -h or --help invoked //
		if(params.containsKey("help") || params.containsKey("h"))
		{
			System.err.println(usage);
			return;
		}

	    // Define the statistics hash and reset counters //


		long numBases = 0;
		long numVariantPositions = 0;
		long numSNPpositions = 0;
		long numIndelPositions = 0;
		long numFailStrandFilter = 0;
		long numVariantsReported = 0;
		long numSNPsReported = 0;
		long numIndelsReported = 0;


		int numParsingExceptions = 0;

		// Parse piped input or user-provided pileup file //

	    try
	    {

	    	// Declare file-parsing variables //

	    	BufferedReader in = VarScan.getInfile(args);
	    	String line;

	    	// If no input, print usage //

	    	if(in == null)
	    	{
	    		System.out.println(usage);
				return;
	    	}

	    	// If input file not ready, give it a few seconds //
	    	int numNaps = 0;

	    	while(!in.ready())
	    	{
	    		try {
			    	Thread.sleep(5000);
			    	numNaps++;

			    	if(numNaps > 100)
			    	{
			    		System.err.println("Input file was not ready after 100 5-second cycles!");
			    		System.exit(10);
			    	}
		    	}
		    	catch(Exception e)
		    	{
		    		System.err.println("Exception while trying to get input" + e.getMessage());
		    		System.exit(1);
		    	}
	    	}

	    	// Proceed if input stream is ready //
			String vcfHeader = "##fileformat=VCFv4.1";

	    	if(in != null && in.ready())
	    	{
	    		// Print a file header //
	    		if(!params.containsKey("no-headers"))
	    		{
	    			if(params.containsKey("output-vcf"))
	    			{
	    				// Output VCF Header //

	    				vcfHeader += "\n" + "##source=VarScan2";
	    				vcfHeader += "\n" + "##INFO=<ID=ADP,Number=1,Type=Integer,Description=\"Average per-sample depth of bases with Phred score >= " + minAvgQual + "\">";
	    				vcfHeader += "\n" + "##INFO=<ID=WT,Number=1,Type=Integer,Description=\"Number of samples called reference (wild-type)\">";
	    				vcfHeader += "\n" + "##INFO=<ID=HET,Number=1,Type=Integer,Description=\"Number of samples called heterozygous-variant\">";
	    				vcfHeader += "\n" + "##INFO=<ID=HOM,Number=1,Type=Integer,Description=\"Number of samples called homozygous-variant\">";
	    				vcfHeader += "\n" + "##INFO=<ID=NC,Number=1,Type=Integer,Description=\"Number of samples not called\">";
	    				vcfHeader += "\n" + "##FILTER=<ID=str10,Description=\"Less than 10% or more than 90% of variant supporting reads on one strand\">";
	    				vcfHeader += "\n" + "##FILTER=<ID=indelError,Description=\"Likely artifact due to indel reads at this position\">";
	    				vcfHeader += "\n" + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">";
	    				vcfHeader += "\n" + "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">";
	    				vcfHeader += "\n" + "##FORMAT=<ID=SDP,Number=1,Type=Integer,Description=\"Raw Read Depth as reported by SAMtools\">";
	    				vcfHeader += "\n" + "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Quality Read Depth of bases with Phred score >= " + minAvgQual + "\">";
	    				vcfHeader += "\n" + "##FORMAT=<ID=RD,Number=1,Type=Integer,Description=\"Depth of reference-supporting bases (reads1)\">";
	    				vcfHeader += "\n" + "##FORMAT=<ID=AD,Number=1,Type=Integer,Description=\"Depth of variant-supporting bases (reads2)\">";
	    				vcfHeader += "\n" + "##FORMAT=<ID=FREQ,Number=1,Type=String,Description=\"Variant allele frequency\">";
	    				vcfHeader += "\n" + "##FORMAT=<ID=PVAL,Number=1,Type=String,Description=\"P-value from Fisher's Exact Test\">";
	    				vcfHeader += "\n" + "##FORMAT=<ID=RBQ,Number=1,Type=Integer,Description=\"Average quality of reference-supporting bases (qual1)\">";
	    				vcfHeader += "\n" + "##FORMAT=<ID=ABQ,Number=1,Type=Integer,Description=\"Average quality of variant-supporting bases (qual2)\">";
	    				vcfHeader += "\n" + "##FORMAT=<ID=RDF,Number=1,Type=Integer,Description=\"Depth of reference-supporting bases on forward strand (reads1plus)\">";
	    				vcfHeader += "\n" + "##FORMAT=<ID=RDR,Number=1,Type=Integer,Description=\"Depth of reference-supporting bases on reverse strand (reads1minus)\">";
	    				vcfHeader += "\n" + "##FORMAT=<ID=ADF,Number=1,Type=Integer,Description=\"Depth of variant-supporting bases on forward strand (reads2plus)\">";
	    				vcfHeader += "\n" + "##FORMAT=<ID=ADR,Number=1,Type=Integer,Description=\"Depth of variant-supporting bases on reverse strand (reads2minus)\">";

	    			}
	    			else
	    			{
	    				// Output VarScan Header //
		    			System.out.println("Chrom\tPosition\tRef\tVar\tCons:Cov:Reads1:Reads2:Freq:P-value\tStrandFilter:R1+:R1-:R2+:R2-:pval\tSamplesRef\tSamplesHet\tSamplesHom\tSamplesNC\tCons:Cov:Reads1:Reads2:Freq:P-value");
	    			}

	    		}



	    		// Parse the infile line by line //

	    		while ((line = in.readLine()) != null)
	    		{
	    			numBases++;//stats.put("numBases", (stats.get("numBases") + 1));

	    			// Output progress line //
	    			if(params.containsKey("verbose") && (numBases % 100000) == 0)
		        		System.err.println(numBases + " positions parsed...");

	    			// Begin try-catch for line parsing //

	    			try
	    			{
	    				String[] lineContents = line.split("\t", -1);

	    				// Output VCF header if desired //

	    				if(numBases == 1 && params.containsKey("output-vcf"))
    	    			{
    	    				vcfHeader += "\n" + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    	    				if(sampleList.length() > 0)
    	    				{
    	    					vcfHeader += "\t" + sampleList;
    	    				}
    	    				else
    	    				{
	    	    				// print the VCF sample header //
	    	    				int sampleCounter = 0;
	    	    				for(int colCounter = 3; colCounter <= (lineContents.length - 3); colCounter += 3)
	    	    				{
	    	    					sampleCounter++;
	    	    					vcfHeader += "\tSample" + sampleCounter;
	    	    				}

    	    				}

    	    				System.out.println(vcfHeader);
    	    			}

	    				// Verify expected pileup format //

	    				if(lineContents.length > 5 && lineContents[0].length() > 0 && lineContents[1].length() > 0 && lineContents[2].length() > 0 && lineContents[3].length() > 0)
	    				{
	    					String refName = lineContents[0];
	    					String position = lineContents[1];
	    					String refBase = lineContents[2].toUpperCase();
	    					String callDepths = "";
	    					String callResults = "";
	    					String vcfResults = "";
	    					HashMap<String, Integer> varAlleles = new HashMap<String, Integer>();
	    					boolean variantFlag = false;
	    					boolean snpFlag = false;
	    					boolean indelFlag = false;
	    					int samplesRef = 0;
	    					int samplesHet = 0;
	    					int samplesHom = 0;
	    					int samplesUncalled = 0;

	    					// Declare variables for cross-sample calling and strand filter //
	    					int allReadDepth = 0;
	    					int allReads1plus = 0;
	    					int allReads1minus = 0;
	    					int allReads2plus = 0;
	    					int allReads2minus = 0;
	    					double strandPvalue = 1.00;
	    					String allReadBases = "";
	    					String allReadQualities = "";

	    					// Call Individual Genotypes for All Samples in Mpileup //

	    					for(int colCounter = 3; colCounter <= (lineContents.length - 3); colCounter += 3)
	    					{
		    					int readDepth = 0;
		    					String readBases = "";
		    					String readQualities = "";
		    					String mapQualities = "";

			    	        	readDepth = Integer.parseInt(lineContents[colCounter]);
			    	        	readBases = lineContents[colCounter + 1];
			    	        	readQualities = lineContents[colCounter + 2];

			    	        	// Append to our long-running total //

			    	        	allReadDepth += readDepth;
			    	        	allReadBases = allReadBases + readBases;
			    	        	allReadQualities = allReadQualities + readQualities;

			    	        	// Determine if this sample's depth meets our minimum //
			    	        	int qualityDepth = 0;
		    	        		qualityDepth = VarScan.qualityDepth(readQualities, minAvgQual);

		    	        		String thisCall = "N" + ":" + qualityDepth + ":-:-:-:-";
		    	        		String thisVCF = "./.:.:" + qualityDepth;

			    	        	if(readDepth >= minCoverage && qualityDepth >= minCoverage)
			    	        	{
			    	        		HashMap<String, String> readCounts = VarScan.getReadCounts(refBase, readBases, readQualities, minAvgQual, mapQualities);
			    	        		String positionCall = VarScan.callPosition(refBase, readCounts, "CNS", minReads2, minVarFreq, minAvgQual, pValueThreshold, minFreqForHom);

			    	        		if(positionCall.length() > 0)
			    	        		{
			    	        			String[] callLines = positionCall.split("\n");

			    	        			// Go thru each line in resulting call list //
			    	        			for(int lineCounter = 0; lineCounter < callLines.length; lineCounter++)
			    	        			{
			    	        				// Determine type of call that was made //
			    	        				String[] callContents = callLines[lineCounter].split("\t");
			    	        				String consBase = callContents[0];
			    	        				int reads1 = Integer.parseInt(callContents[1]);
			    	        				int reads2 = Integer.parseInt(callContents[2]);
			    	        				String varFreq = callContents[3];
			    	        				int strands1 = Integer.parseInt(callContents[4]);
			    	        				int strands2 = Integer.parseInt(callContents[5]);
			    	        				int qual1 = Integer.parseInt(callContents[6]);
			    	        				int qual2 = Integer.parseInt(callContents[7]);
			    	        				double pValue = Double.parseDouble(callContents[8]);
			    	        				int reads1plus = Integer.parseInt(callContents[11]);
			    	        				int reads1minus = Integer.parseInt(callContents[12]);
			    	        				int reads2plus = Integer.parseInt(callContents[13]);
			    	        				int reads2minus = Integer.parseInt(callContents[14]);
			    	        				String varAllele = "";

			    	        				double logP = 0;
			    	        				try {
				    	        				logP = 0 - (10 * java.lang.Math.log10(pValue));
				    	        				if(logP > 255)
				    	        					logP = 255;
			    	        				}
			    	        				catch(Exception e)
			    	        				{
			    	        					// Stick with default logP value
			    	        				}


			    	        				// Capture the variant allele if there is one //

			    	        				if(!consBase.equals(refBase) && !consBase.equals("N") && callContents.length > 15)
			    	        				{
			    	        					varAllele = callContents[15];

			    	        					// Determine how many variant alleles have been seen //

			    	        					int varAlleleNumber = 0;

			    	        					// Determine if we've seen the variant and what its number is ##

				    	        				if(varAlleles.containsKey(varAllele))
				    	        				{
				    	        					varAlleleNumber = varAlleles.get(varAllele);
				    	        				}
				    	        				else
				    	        				{
				    	        					// IF no variants yet seen, this is variant allele 1 //
				    	        					varAlleleNumber = varAlleles.size() + 1;
				    	        					varAlleles.put(varAllele, varAlleleNumber);
				    	        				}

			    	        					if(VarScan.isHomozygous(consBase))
			    	        					{
			    	        						samplesHom++;
			    	        						thisVCF = varAlleleNumber + "/" + varAlleleNumber;
			    	        					}
			    	        					else
			    	        					{
			    	        						samplesHet++;
			    	        						thisVCF = "0" + "/" + varAlleleNumber;
			    	        					}

			    	        					thisVCF += ":" + (int) logP + ":" + readDepth + ":" + qualityDepth;
		    	        						thisVCF += ":" + reads1 + ":" + reads2 + ":" + varFreq + ":" + pvalueFormat.format(pValue);
		    	        						thisVCF += ":" + qual1 + ":" + qual2;
		    	        						thisVCF += ":" + reads1plus + ":" + reads1minus + ":" + reads2plus + ":" + reads2minus;
			    	        				}
			    	        				else if(consBase.equals(refBase))
			    	        				{
			    	        					// A reference call - recalculate p-value against a possible het //
			    	        					int expReads1 = (reads1 + reads2) / 2;
			    	        					int expReads2 = (reads1 + reads2) - expReads1;
			    	        					double newPvalue = VarScan.getSignificance(reads1, reads2, expReads1, expReads2);
			    	        					double newLogP = 0;
				    	        				try {
					    	        				newLogP = 0 - (10 * java.lang.Math.log10(newPvalue));
				    	        				}
				    	        				catch(Exception e)
				    	        				{
				    	        					// Stick with default logP value
				    	        				}
			    	        					thisVCF = "0" + "/" + "0";
			    	        					thisVCF += ":" + (int) newLogP + ":" + readDepth + ":" + qualityDepth;
		    	        						thisVCF += ":" + reads1 + ":" + reads2 + ":" + varFreq + ":" + pvalueFormat.format(pValue);
		    	        						thisVCF += ":" + qual1 + ":" + qual2;
		    	        						thisVCF += ":" + reads1plus + ":" + reads1minus + ":" + reads2plus + ":" + reads2minus;

			    	        				}


			    	        				thisCall = consBase + ":" + qualityDepth + ":" + reads1 + ":" + reads2 + ":" + varFreq;
			    	        				thisCall += ":" + pvalueFormat.format(pValue);

			    	        				if(!consBase.equals(refBase) && !consBase.equals("N"))
			    	        				{
			    	        					variantFlag = true;

			    	        					// Flag what type of variant was observed //
			    	        					if(consBase.length() > 1)
			    	        						indelFlag = true;
			    	        					else
			    	        						snpFlag = true;

			    	        					// Save reads1plus and reads1minus //

			    	        					allReads1plus += reads1plus;
			    	        					allReads1minus += reads1minus;
			    	        					allReads2plus += reads2plus;
			    	        					allReads2minus += reads2minus;


			    	        				}
			    	        				else
			    	        				{
			    	        					samplesRef++;
			    	        				}
			    	        			}

			    	        		}
			    	        		else
			    	        		{
			    	        			samplesUncalled++;
			    	        		}


			    	        	}
			    	        	else
			    	        	{
			    	        		samplesUncalled++;
			    	        	}


			    	        	// Add this depth to the list //

			    	        	if(callDepths.length() > 0)
			    	        		callDepths += " ";

			    	        	callDepths += readDepth;

		    	        		// Add this call to the list //
		    	        		if(callResults.length() > 0)
		    	        			callResults = callResults + " ";

		    	        		callResults = callResults + thisCall;

		    	        		// Add this to the sample VCF string //

		    	        		if(vcfResults.length() > 0)
		    	        			vcfResults = vcfResults + "\t";

		    	        		vcfResults = vcfResults + thisVCF;
	    					}


	    					// Call the cross-sample pileup //

		    	        	int qualityDepth = 0;
	    	        		qualityDepth = VarScan.qualityDepth(allReadQualities, minAvgQual);
		    	        	String allMapQualities = "";
		    	        	String allConsensusCall = "N:" + qualityDepth + ":-:-:-:-";



		    	        	if(params.containsKey("output-vcf"))
		    	        	{
		    	        		// Skip this if we're outputting VCF //

		    	        	}
		    	        	else if(allReadDepth >= minCoverage && qualityDepth >= minCoverage)
		    	        	{
		    	        		HashMap<String, String> readCounts = VarScan.getReadCounts(refBase, allReadBases, allReadQualities, minAvgQual, allMapQualities);
		    	        		String positionCall = VarScan.callPosition(refBase, readCounts, "CNS", minReads2, minVarFreq, minAvgQual, pValueThreshold, minFreqForHom);

		    	        		if(positionCall.length() > 0)
		    	        		{
		    	        			String[] callLines = positionCall.split("\n");

		    	        			// Go thru each line in resulting call list //
		    	        			for(int lineCounter = 0; lineCounter < callLines.length; lineCounter++)
		    	        			{
		    	        				// Determine type of call that was made //
		    	        				String[] callContents = callLines[lineCounter].split("\t");
		    	        				String consBase = callContents[0];
		    	        				int reads1 = Integer.parseInt(callContents[1]);
		    	        				int reads2 = Integer.parseInt(callContents[2]);
		    	        				String varFreq = callContents[3];
		    	        				double pValue = Double.parseDouble(callContents[8]);
		    	        				String varAllele = "";

		    	        				// Capture the variant allele if there is one //

		    	        				if(!consBase.equals(refBase) && callContents.length > 15)
		    	        				{
		    	        					varAllele = callContents[15];
		    	        					if(varAlleles.containsKey(varAllele))
			    	        				{
//			    	        					varAlleles.put(varAllele, (varAlleles.get(varAllele) + 1));
			    	        				}
			    	        				else
			    	        				{
			    	        					// IF no variants yet seen, this is variant allele 1 //
			    	        					int varAlleleNumber = varAlleles.size() + 1;
			    	        					varAlleles.put(varAllele, varAlleleNumber);
			    	        				}

		    	        				}


		    	        				allConsensusCall = consBase + ":" + qualityDepth + ":" + reads1 + ":" + reads2 + ":" + varFreq;
		    	        				allConsensusCall += ":" + pvalueFormat.format(pValue);

		    	        				if(!consBase.equals(refBase) && !consBase.equals("N"))
		    	        				{
		    	        					variantFlag = true;

		    	        					// Flag what type of variant was observed //
		    	        					if(consBase.length() > 1)
		    	        						indelFlag = true;
		    	        					else
		    	        						snpFlag = true;

		    	        				}
		    	        			}

		    	        		}
		    	        		else
		    	        		{
		    	        			// NO call made from all-sample pileup //
		    	        		}


		    	        	}
		    	        	else
		    	        	{
		    	        		// All-sample pileup failed to meet min depth //
		    	        	}


	    					// Get All Variant alleles observed //

	    					String varBases = "";
	    					// First, obtain their unique keys which are in alphanumeric order //
	    					String[] sortedKeys = (String[]) varAlleles.keySet().toArray(new String[0]);

	    					// Create an empty array to put these into sorted order //
	    					String[] alleleKeys = new String[sortedKeys.length];

	    					// Put alleles into this array in their order of occurrence in VCF line //
	    					for(String allele : sortedKeys)
	    					{
	    						int arrayIndex = varAlleles.get(allele) - 1;
    							alleleKeys[arrayIndex] = allele;
	    					}

	    					// Export all variant alleles into a comma-separated string//
	    					// This is what's provided in native output, or converted to VCF format //
	    					for(String allele : alleleKeys)
	    					{
	    						if(varBases.length() > 0)
	    							varBases += ",";

	    						varBases += allele;
	    					}

	    					// It's possible that we see no variant here, so we need the proper empty character //
	    					if(varBases.length() == 0)
	    						varBases = ".";

	    					// Count whether there was a variant //
    						if(variantFlag)
    							numVariantPositions++;
    						if(snpFlag)
    							numSNPpositions++;
    						if(indelFlag)
    							numIndelPositions++;

	    					// Determine strand filter status if it's turned on  //
	    					String strandFilterStatus = "Pass:" + allReads1plus + ":" + allReads1minus + ":" + allReads2plus + ":" + allReads2minus + ":" + pvalueFormat.format(strandPvalue);
	    					boolean failedStrandFilter = false;

    						if(strandFilter && variantFlag && (allReads1plus > 0 || allReads1minus > 0 || allReads2plus > 0 || allReads2plus > 0))
	    					{
    	    					double refStrandPlus = 0.50;
    	    					double varStrandPlus = 0.50;

    	    					// Calculate strandedness for variant allele //

    	    					if((allReads2plus + allReads2minus) > 0)
    	    						varStrandPlus = (double) allReads2plus / (double) (allReads2plus + allReads2minus);

    	    					// To save time, only calculate p-value if var strandedness is biased //

    	    					if(varStrandPlus < 0.10 || varStrandPlus > 0.90)
    	    					{
    	    						// Calculate strandedness for reference allele if we have 2+ reads //

        	    					if((allReads1plus + allReads1minus) > 1)
        	    					{
        	    						refStrandPlus = (double) allReads1plus / (double) (allReads1plus + allReads1minus);
        	    						strandPvalue = VarScan.getSignificance(allReads1plus, allReads1minus, allReads2plus, allReads2minus);
        	    					}
        	    					// Otherwise, only homozygous-variant reads seen, so compare to a 50/50 distribution //
    	    						else
    	    						{
    	    							// Compare to expected 50/50 distribution //
    	    							int testReads1plus = (int) (allReads2plus + allReads2minus) / 2;
    	    							int testReads1minus = (allReads2plus + allReads2minus) - testReads1plus;
    	    							strandPvalue = VarScan.getSignificance(testReads1plus, testReads1minus, allReads2plus, allReads2minus);
    	    						}

        	    					strandFilterStatus = "Pass:" + varStrandPlus + ":" + allReads1plus + ":" + allReads1minus + ":" + allReads2plus + ":" + allReads2minus + ":" + pvalueFormat.format(strandPvalue);

    	    						// If ref allele had good strandedness, and var allele did not, this may be a failure //
    	    						if(refStrandPlus >= 0.10 && refStrandPlus <= 0.90 && !(varStrandPlus >= 0.10 && varStrandPlus <= 0.90))
    	    						{
    	    							if(strandPvalue < strandPvalueThreshold)
    	    	    					{
    	    	    						strandFilterStatus = "Fail:" + allReads1plus + ":" + allReads1minus + ":" + allReads2plus + ":" + allReads2minus + ":" + pvalueFormat.format(strandPvalue);
    	    	    						numFailStrandFilter++;
    	    	    						failedStrandFilter = true;
    	    	    					}
    	    						}
    	    					}

	    					}

    						// BEGIN BUILDING OUTPUT LINE //

	    					String outLine = refName + "\t" + position + "\t";

	    					if(params.containsKey("output-vcf"))
	    					{
	    						// Calculate average sample depth //
	    						int avgQualityDepth = 0;
	    						if((samplesRef + samplesHet + samplesHom + samplesUncalled) > 0)
	    							avgQualityDepth = qualityDepth / (samplesRef + samplesHet + samplesHom + samplesUncalled);

	    						String refColumn = "";
	    						String varColumn = "";

	    						// Handle complex positions with multiple alleles including at least one indel //

	    						if(varBases.contains(",") && (varBases.contains("-") || varBases.contains("+")))
	    						{
	    							// Multi-allele indel //
	    							int maxDelSize = 0;
	    							String maxDelBases = "";
	    							// Go through each varAllele to find longest deletion //
	    							String[] varBaseContents = varBases.split(",");
	    							for(String varAllele : varBaseContents)
	    							{
	    								if(varAllele.startsWith("-"))
	    								{
	    									varAllele = varAllele.replace("-", "");
	    									if(varAllele.length() > maxDelSize)
	    									{
	    										maxDelBases = varAllele;
	    										maxDelSize = varAllele.length();
	    									}
	    								}
	    							}

	    							// Set refBase to maximum del //
	    							refColumn = refBase + maxDelBases;

	    							// Establish each allele in var Column //
	    							varColumn = "";

	    							for(String varAllele : varBaseContents)
	    							{
    									if(varColumn.length() > 0)
    										varColumn = varColumn + ",";

	    								if(varAllele.startsWith("-"))
	    								{
	    									varAllele = varAllele.replace("-", "");

	    									// For the smaller deletion, determine ref bases to add //
	    									if(varAllele.length() < maxDelSize)
	    									{
	    										String varEntry = maxDelBases.replace(varAllele, "");
	    										varColumn = varColumn + refBase + varEntry;
	    									}
	    									else
	    									{
	    										varColumn = varColumn + refBase;
	    									}
	    								}
	    								else if(varAllele.startsWith("+"))
	    								{
	    									varAllele = varAllele.replace("+", "");
	    									String varEntry = refBase + varAllele + maxDelBases;
	    									varColumn = varColumn + varEntry;
	    								}
	    								else
	    								{
	    									String varEntry = varAllele + maxDelBases;
	    									varColumn = varColumn + varEntry;
	    								}
	    							}


	    						}

	    						else if(varBases.startsWith("+"))
		    					{
		    						// INSERTION //
		    						// Ref = ref base; Var = ref base followed by inserted bases //
	    							refColumn = refBase;
		    						varColumn = refBase + varBases.replace("+", "");
		    					}
		    					else if(varBases.startsWith("-"))
		    					{
		    						// DELETION //
		    						// Ref = ref base followed by deleted bases; var = ref base //
		    						refColumn = refBase + varBases.replace("-", "");
		    						varColumn = refBase;
		    					}
		    					else
		    					{
		    						refColumn = refBase;
		    						varColumn = varBases;
		    					}

	    						// Ensure that varColumn does not contain any +/- //
	    						varColumn = varColumn.replace("+", "");
	    						varColumn = varColumn.replace("-", "");


	    						outLine += "." + "\t" + refColumn + "\t" + varColumn + "\t.\t";

	    						if(strandFilterStatus.contains("Pass"))
	    							outLine += "PASS\t";
	    						else
	    							outLine += "str10\t";
	    						outLine += "ADP=" + avgQualityDepth + ";WT=" + samplesRef + ";HET=" + samplesHet + ";HOM=" + samplesHom + ";NC=" + samplesUncalled;
	    						outLine += "\t" + "GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR" + "\t";
	    						outLine += vcfResults;
	    					}
	    					else
	    					{
	    						outLine += refBase + "\t" + varBases + "\t";
		    					outLine += allConsensusCall + "\t" + strandFilterStatus + "\t";
		    					outLine += samplesRef + "\t" + samplesHet + "\t" + samplesHom + "\t" + samplesUncalled + "\t";
		    					outLine += callResults;
	    					}



    						// If there was a variant, but strand-filter failed, and only reporting variants:
	    					boolean reportFlag = false;

    						if(variantFlag && strandFilter && failedStrandFilter)
    						{
    							// Do not print a variant that failed strand-filter unless in CNS mode //
    							if(!variantsOnly && !snpsOnly && !indelsOnly)
    								reportFlag = true;
    						}
    						else if((variantsOnly || snpsOnly || indelsOnly) && !variantFlag)
    						{
    							// Do not print if reporting variants, but no variant was seen //
    						}
    						else if(!variantsOnly && !snpsOnly && !indelsOnly)
	    					{
    							// Print consensus if in consensus calling mode //
    							reportFlag = true;
	    					}
    						else if(variantFlag && variantsOnly)
    						{
    							// Print any variant if variants flag set //
    							reportFlag = true;
    						}
    						else if(snpFlag && snpsOnly)
    						{
    							// Print SNP variant if SNPs-only flag set //
    							reportFlag = true;
    						}
    						else if(indelFlag && indelsOnly)
    						{
    							// Print indel variant if indels-only flag set //
    							reportFlag = true;
    						}
    						else
    						{
    							// Don't report a consensus call if limited to variant reporting //
    						}

    						if(reportFlag)
    						{
	    						System.out.println(outLine);

	    						if(variantFlag)
	    							numVariantsReported++;
	    						if(snpFlag)
	    							numSNPsReported++;
	    						if(indelFlag)
	    							numIndelsReported++;
    						}

	    				}
	    				else
	    				{
	    					if(lineContents.length >= 4 && lineContents[3].equals("0"))
	    					{
	    						// A pileup line with 0x coverage, so ignore
	    					}
	    					else
	    					{
		    					System.err.println("Warning: Line ignored: Invalid format for pileup at line " + numBases + "\n" + line + "\n");
		    					return;
	    					}

	    				}
	    			}
	    			catch(Exception e)
	    		    {
	    		    	System.err.println("Parsing Exception on line:\n" + line + "\n" + e.getLocalizedMessage());
	    				numParsingExceptions++;
	    				if(numParsingExceptions >= 5)
	    				{
	    					System.err.println("Too many parsing exceptions encountered; exiting");
	    					return;
	    				}
	    		    	return;
	    		    }


	    		}

				in.close();

				System.err.println(numBases + " bases in pileup file");
				System.err.println(numVariantPositions + " variant positions (" + numSNPpositions + " SNP, " + numIndelPositions + " indel)");
				System.err.println(numFailStrandFilter + " were failed by the strand-filter");
				System.err.println(numVariantsReported + " variant positions reported (" + numSNPsReported + " SNP, " + numIndelsReported + " indel)");
	    	}
	    	// Insufficient input was provided, so print usage //
	    	else
	    	{
				 System.err.println("Please provide an input file!\n" + usage);
				 System.exit(10);
	    	}
	    }
	    catch(Exception e)
	    {
	    	System.err.println("Exception: " + e.getLocalizedMessage());
	    	e.printStackTrace(System.err);
	    	System.exit(11);
	    }
	}
}
