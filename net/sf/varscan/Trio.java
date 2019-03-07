/**
 * @(#)Trio.java
 *
 * Copyright (c) 2009-2013 Daniel C. Koboldt and Washington University in St. Louis
 *
 * COPYRIGHT
 */

package net.sf.varscan;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.PrintStream;
import java.text.DecimalFormat;
import java.util.*;
import java.lang.Math;

/**
 * A class for calling variants in a mother-father-child trio
 *
 * @version	2.3
 *
 * @author Daniel C. Koboldt <dkoboldt@genome.wustl.edu>
 *
 */
public class Trio {
	public Trio(String[] args, String callType)
	{
		// Define the usage message //
		String usage = "USAGE: java -jar VarScan.jar trio [mpileup file] [output-basename] OPTIONS\n" +
		"\tmpileup file - The SAMtools mpileup file for father, mother, child in that order\n" +
		"\n" +
		"\tOPTIONS:\n" +
		"\t--output-name\tAn output base name for VCF files of results. Required for piped input\n" +
		"\t--min-coverage\tMinimum read depth at a position to make a call [20]\n" +
		"\t--min-reads2\tMinimum supporting reads at a position to call variants [2]\n" +
		"\t--min-avg-qual\tMinimum base quality at a position to count a read [15]\n" +
		"\t--min-var-freq\tMinimum variant allele frequency threshold [0.20]\n" +
		"\t--min-freq-for-hom\tMinimum frequency to call homozygote [0.75]\n" +
		"\t--p-value\tDefault p-value threshold for calling variants [0.05]\n" +
		"\t--adj-var-freq\tAdjusted minimum VAF when recalling at variant site [0.05]\n" +
		"\t--adj-p-value\tAdjusted p-value when recalling at variant site [0.10]\n" +
		"\t--vcf-sample-list\tFor VCF output, a list of sample names in order, one per line\n" +
		"\t--variants\tReport only variant (SNP/indel) positions [0]";

		// Set parameter defaults //

		HashMap<String, String> params = VarScan.getParams(args);

		// Establish output file names //
		String outputName = "output";
		String outputSnp = "";
		String outputIndel = "";

		if(args.length >= 3 && !args[2].startsWith("-"))
		{
			outputName = args[2];
			outputSnp = outputName + ".snp.vcf";
			outputIndel = outputName + ".indel.vcf";
		}

		// Set up formatting for p-values //
		DecimalFormat pvalueFormat = new DecimalFormat("0.####E0");

		// Force VCF output //
		params.put("output-vcf", "1");

		// Set parameter defaults //

		int minCoverage = 20;
		int minReads2 = 4;
		int minAvgQual = 15;
		double minVarFreq = 0.20;
		double minFreqForHom = 0.75;
		double pValueThreshold = 0.01;
		double strandPvalueThreshold = 0.01;
		double adjustedMinVarFreq = 0.05;
		double adjustedpValueThreshold = 0.10;
		int adjustedMinReads2 = 2;
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

			if(params.containsKey("adj-min-reads2"))
				 adjustedMinReads2 = Integer.parseInt(params.get("adj-min-reads2"));

			if(params.containsKey("min-var-freq"))
				 minVarFreq = Double.parseDouble(params.get("min-var-freq"));

			if(params.containsKey("adj-var-freq"))
				 adjustedMinVarFreq = Double.parseDouble(params.get("adj-var-freq"));

			if(params.containsKey("p-value"))
				 pValueThreshold = Double.parseDouble(params.get("p-value"));

			if(params.containsKey("adj-p-value"))
				 adjustedpValueThreshold = Double.parseDouble(params.get("adj-p-value"));

			if(params.containsKey("min-freq-for-hom"))
				 minFreqForHom = Double.parseDouble(params.get("min-freq-for-hom"));

			if(params.containsKey("min-avg-qual"))
				 minAvgQual = Integer.parseInt(params.get("min-avg-qual"));

			if(params.containsKey("output-name"))
			{
				outputName = params.get("output-name");
				outputSnp = outputName + ".snp.vcf";
				outputIndel = outputName + ".indel.vcf";

			}

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

			if(params.containsKey("p-value"))
				 pValueThreshold = Double.parseDouble(params.get("p-value"));
			else
				System.err.println("Warning: No p-value threshold provided, so p-values will not be calculated");


			// Check for correct input //

			if(outputSnp.length() == 0 || outputIndel.length() == 0)
			{
				System.err.println("Please provide an output basename or SNP/indel output files!");
				 System.err.println(usage);
				 System.exit(1);
			}

			System.err.println("SNPs will be output to " + outputSnp);
			System.err.println("Indels will be output to " + outputIndel);
			System.err.println("Min coverage:\t" + minCoverage);
			 System.err.println("Min reads2:\t" + minReads2);
			 System.err.println("Min var freq:\t" + minVarFreq);
			 System.err.println("Min avg qual:\t" + minAvgQual);
			 System.err.println("P-value thresh:\t" + pValueThreshold);
			 System.err.println("Adj. min reads2:\t" + adjustedMinReads2);
			 System.err.println("Adj. var freq:\t" + adjustedMinVarFreq);
			 System.err.println("Adj. p-value:\t" + adjustedpValueThreshold);
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
		long numBasesCovered = 0;
		long numVariantPositions = 0;
		long numSNPpositions = 0;
		long numIndelPositions = 0;
		long numFailStrandFilter = 0;
		long numFailMendelFilter = 0;
		long numVariantsReported = 0;
		long numVariantsReportedDeNovo = 0;
		long numSNPsReported = 0;
		long numSNPsReportedDeNovo = 0;
		long numIndelsReported = 0;
		long numIndelsReportedDeNovo = 0;

		int numParsingExceptions = 0;

		HashMap<String, Integer> stats = new HashMap<String, Integer>();

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
				// Declare output file //
		 	 	PrintStream outSnp = null; // declare a print stream object for SNPs
		 	 	PrintStream outIndel = null; // declare a print stream object for Indels

		 		outSnp = new PrintStream( new FileOutputStream(outputSnp) );
		 		outIndel = new PrintStream( new FileOutputStream(outputIndel) );

	    		// Print a file header //
	    		if(!params.containsKey("no-headers"))
	    		{
	    			if(params.containsKey("output-vcf"))
	    			{
	    				// Output VCF Header //

	    				vcfHeader += "\n" + "##source=VarScan2";
	    				vcfHeader += "\n" + "##INFO=<ID=ADP,Number=1,Type=Integer,Description=\"Average per-sample depth of bases with Phred score >= " + minAvgQual + "\">";
	    				vcfHeader += "\n" + "##INFO=<ID=STATUS,Number=1,Type=String,Description=\"Variant status in trio (1=untransmitted, 2=transmitted, 3=denovo, 4=MIE)\">";
	    				vcfHeader += "\n" + "##INFO=<ID=DENOVO,Number=0,Type=Flag,Description=\"Indicates apparent de novo mutations unique to the child\">";
	    				vcfHeader += "\n" + "##FILTER=<ID=str10,Description=\"Less than 10% or more than 90% of variant supporting reads on one strand\">";
	    				vcfHeader += "\n" + "##FILTER=<ID=indelError,Description=\"Likely artifact due to indel reads at this position\">";
	    				vcfHeader += "\n" + "##FILTER=<ID=mendelError,Description=\"Apparent Mendelian inheritance error (MIE) in trio\">";
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
		    			System.out.println("Chrom\tPosition\tRef\tVar\tStrandFilter:R1+:R1-:R2+:R2-:pval\tFather:Cov:Reads1:Reads2:Freq:P-value\tMother:Cov:Reads1:Reads2:Freq:P-value\tChild:Cov:Reads1:Reads2:Freq:P-value");
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
	    				String[] lineContents = line.split("\t");

	    				// Verify expected pileup format //

	    				if(lineContents.length > 5 && lineContents[0].length() > 0 && lineContents[1].length() > 0 && lineContents[2].length() > 0 && lineContents[3].length() > 0)
	    				{
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
									vcfHeader += "\tFather\tMother\tChild";
	    	    				}

	    						// Output VCF Header //
	    						outSnp.println(vcfHeader);
	    						outIndel.println(vcfHeader);

	    	    			}


	    					String refName = lineContents[0];
	    					String position = lineContents[1];
	    					String refBase = lineContents[2].toUpperCase();
	    					HashMap<String, Integer> varAlleles = new HashMap<String, Integer>();
	    					boolean variantFlag = false;

	    					// Declare variables for cross-sample calling and strand filter //
	    					double strandPvalue = 1.00;
	    					String strandFilterStatus = "";

	    					if(lineContents.length > 12)
	    					{
	    						if(numBases == 1)
	    							System.err.println("Warning: More than 3 samples in pileup; but only first 3 will be used and they should be father, mother child");
	    					}

	    					// Get Father Call //
	    					int offset = 3;
	    					int fatherDepth = Integer.parseInt(lineContents[offset]);
		    	        	String fatherBases = lineContents[offset + 1];
		    	        	String fatherQualities = lineContents[offset + 2];
		    	        	int fatherQualityDepth = VarScan.qualityDepth(fatherQualities, minAvgQual);

	    	        		// Get Mother Call //
	    					offset = 6;
	    					int motherDepth = Integer.parseInt(lineContents[offset]);
		    	        	String motherBases = lineContents[offset + 1];
		    	        	String motherQualities = lineContents[offset + 2];
		    	        	int motherQualityDepth = VarScan.qualityDepth(motherQualities, minAvgQual);

	    	        		// Get Child Call //
	    					offset = 9;
	    					int childDepth = Integer.parseInt(lineContents[offset]);
		    	        	String childBases = lineContents[offset + 1];
		    	        	String childQualities = lineContents[offset + 2];
		    	        	int childQualityDepth = VarScan.qualityDepth(childQualities, minAvgQual);

		    	        	if(fatherQualityDepth >= minCoverage && motherQualityDepth >= minCoverage && childQualityDepth >= minCoverage)
	    	        		{
	    	        			numBasesCovered++;

	    						// Perform strand filter test //
	    						String allBases = fatherBases + motherBases + childBases;
			    	        	String allQualities = fatherQualities + motherQualities + childQualities;
			    	        	HashMap<String, String> allCounts = VarScan.getReadCounts(refBase, allBases, allQualities, minAvgQual, "");
			    	        	String positionCall = VarScan.callPosition(refBase, allCounts, "CNS", minReads2, 0.01, minAvgQual, 0.95, minFreqForHom);
			    	        	String[] callContents = positionCall.split("\t");
			    	        	if(callContents.length >= 15)
			    	        	{
			        				int reads1plus = Integer.parseInt(callContents[11]);
			        				int reads1minus = Integer.parseInt(callContents[12]);
			        				int reads2plus = Integer.parseInt(callContents[13]);
			        				int reads2minus = Integer.parseInt(callContents[14]);
			        				strandFilterStatus = VarScan.strandFilter(reads1plus, reads1minus, reads2plus, reads2minus, strandPvalueThreshold);
			    	        	}
//		        				System.err.println(strandFilterStatus);


	    	        			HashMap<String, String> fatherCounts = VarScan.getReadCounts(refBase, fatherBases, fatherQualities, minAvgQual, "");
		    	        		HashMap<String, String> motherCounts = VarScan.getReadCounts(refBase, motherBases, motherQualities, minAvgQual, "");
		    	        		HashMap<String, String> childCounts = VarScan.getReadCounts(refBase, childBases, childQualities, minAvgQual, "");

		    	        		// Prepare Strings for Results //
		    	        		String fatherCall = "";
		    	        		String motherCall = "";
		    	        		String childCall = "";
		    	        		String trioStatus = "";

		    	        		// Try trio calling //
		    	        		String trioCall = callTrio(refBase, fatherCounts, motherCounts, childCounts, minReads2, minVarFreq, minAvgQual, pValueThreshold, minFreqForHom);
		    	        		String[] trioCallContents = trioCall.split("\t");

		    	        		if(trioCallContents.length >= 4)
		    	        		{
		    	        			fatherCall = trioCallContents[0];
		    	        			motherCall = trioCallContents[1];
		    	        			childCall = trioCallContents[2];

		    	        			trioStatus = trioCallContents[3];
	    	        				boolean recallTrio = false;

		    	        			// Consider re-calling de novo mutations, MIEs, and untransmitted //
		    	        			if(trioStatus.equals("DeNovo") || trioStatus.contains("MIE") || trioStatus.equals("Untransmitted"))
		    	        			{
		    	        				// Parse out the variant allele from each sample //
		    	        				String[] fatherContents = fatherCall.split(":");
		    	        				String[] motherContents = motherCall.split(":");
		    	        				String[] childContents = childCall.split(":");

		    	        				String fatherAllele = refBase;
		    	        				String motherAllele = refBase;
		    	        				String childAllele = refBase;

		    	        				if(fatherContents.length >= 16)
		    	        					fatherAllele = fatherContents[15];

		    	        				if(motherContents.length >= 16)
		    	        					motherAllele = motherContents[15];

		    	        				if(childContents.length >= 16)
		    	        					childAllele = childContents[15];

		    	        				// Evaluate if we should re-call the trio with reduced thresholds //

		    	        				if(trioStatus.equals("Untransmitted"))
		    	        				{
		    	        					// Re-call if child was Reference but has evidence of variant //
		    	        					if(!childAllele.equals(refBase) && (childAllele.equals(motherAllele) || childAllele.equals(fatherAllele)))
		    	        						recallTrio = true;
		    	        				}
		    	        				else if(trioStatus.equals("DeNovo"))
		    	        				{
		    	        					// Recall if child had de novo but either parent has evidence //
		    	        					//if(fatherAllele.equals(childAllele) || motherAllele.equals(childAllele))
		    	        						recallTrio = true;
		    	        				}
		    	        				else if(trioStatus.contains("MIE"))
		    	        				{
		    	        					// Recall if there was an apparent inheritance error //
		    	        					recallTrio = true;
		    	        				}

		    	        				if(recallTrio)
		    	        				{
		    	        					// Adjust values and recall trio //
//		    	        					double adjustedMinVarFreq = minVarFreq / 2.00;
//			    	        				double adjustedpValueThreshold = 0.20;

			    	        				trioCall = callTrio(refBase, fatherCounts, motherCounts, childCounts, adjustedMinReads2, adjustedMinVarFreq, minAvgQual, adjustedpValueThreshold, minFreqForHom);
			    	        				trioCallContents = trioCall.split("\t");

			    	        				// Determine if something changed //
			    	        				if(!trioStatus.equals(trioCallContents[3]))
			    	        				{
			    	        					String change = "initially " + trioStatus + " were re-called " + trioCallContents[3];
					    						if(!stats.containsKey(change))
					    						{
//					    							System.err.println("CHANGED FROM " + trioStatus + "\t" + trioCall);
					    							stats.put(change, 1);
					    						}
					    						else
					    						{
					    							stats.put(change, (stats.get(change) + 1));
					    						}
			    	        				}

			    	        				trioStatus = trioCallContents[3];

		    	        				}	// Otherwise don't re-call //

		    	        			}
		    	        			else
		    	        			{
		    	        				// Must have been Reference or Germline //
		    	        			}

		    	        			String variantType = "SNP";

		    	        			fatherCall = trioCallContents[0];
		    	        			motherCall = trioCallContents[1];
		    	        			childCall = trioCallContents[2];

	    	        				// Parse out the variant allele from each sample //
	    	        				String[] fatherContents = fatherCall.split(":");
	    	        				String[] motherContents = motherCall.split(":");
	    	        				String[] childContents = childCall.split(":");

	    	        				String fatherAllele = refBase;
	    	        				String motherAllele = refBase;
	    	        				String childAllele = refBase;

	    	        				// BUILD FATHER VCF //

	    	        				String fatherVCF = "./.:.:" + fatherQualityDepth;

	    	        				if(fatherContents.length >= 15)
	    	        				{
	    	        					if(fatherContents.length >= 16)
	    	        						fatherAllele = fatherContents[15];
	    	        					String consBase = fatherContents[0];
		    	        				int reads1 = Integer.parseInt(fatherContents[1]);
		    	        				int reads2 = Integer.parseInt(fatherContents[2]);
		    	        				String varFreq = fatherContents[3];
		    	        				int qual1 = Integer.parseInt(fatherContents[6]);
		    	        				int qual2 = Integer.parseInt(fatherContents[7]);
		    	        				double pValue = Double.parseDouble(fatherContents[8]);
		    	        				int reads1plus = Integer.parseInt(fatherContents[11]);
		    	        				int reads1minus = Integer.parseInt(fatherContents[12]);
		    	        				int reads2plus = Integer.parseInt(fatherContents[13]);
		    	        				int reads2minus = Integer.parseInt(fatherContents[14]);

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

		    	        				// Father is wildtype //
		    	        				if(consBase.equals(refBase))
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
		    	        					fatherVCF = "0" + "/" + "0";
		    	        					fatherVCF += ":" + (int) newLogP + ":" + fatherDepth + ":" + fatherQualityDepth;
	    	        						fatherVCF += ":" + reads1 + ":" + reads2 + ":" + varFreq + ":" + pvalueFormat.format(pValue);
	    	        						fatherVCF += ":" + qual1 + ":" + qual2;
	    	        						fatherVCF += ":" + reads1plus + ":" + reads1minus + ":" + reads2plus + ":" + reads2minus;

		    	        				}
		    	        				// Father is variant //
		    	        				else if(fatherAllele.length() > 0 && !fatherAllele.equals("N") && !fatherAllele.equals("."))
		    	        				{
		    	        					// Determine how many variant alleles have been seen //

		    	        					int varAlleleNumber = 0;

		    	        					// Determine if we've seen the variant and what its number is ##

			    	        				if(varAlleles.containsKey(fatherAllele))
			    	        				{
			    	        					varAlleleNumber = varAlleles.get(fatherAllele);
			    	        				}
			    	        				else
			    	        				{
			    	        					// IF no variants yet seen, this is variant allele 1 //
			    	        					varAlleleNumber = varAlleles.size() + 1;
			    	        					varAlleles.put(fatherAllele, varAlleleNumber);
			    	        				}

			    	        				if(fatherContents.length >= 1)
			    	        				{
			    	        					if(VarScan.isHomozygous(consBase))
			    	        					{
			    	        						fatherVCF = varAlleleNumber + "/" + varAlleleNumber;
			    	        					}
			    	        					else
			    	        					{
			    	        						fatherVCF = "0" + "/" + varAlleleNumber;
			    	        					}

			    	        					fatherVCF += ":" + (int) logP + ":" + fatherDepth + ":" + fatherQualityDepth;
		    	        						fatherVCF += ":" + reads1 + ":" + reads2 + ":" + varFreq + ":" + pvalueFormat.format(pValue);
		    	        						fatherVCF += ":" + qual1 + ":" + qual2;
		    	        						fatherVCF += ":" + reads1plus + ":" + reads1minus + ":" + reads2plus + ":" + reads2minus;
			    	        				}

		    	        					if(fatherAllele.length() > 1)
		    	        						variantType = "Indel";
		    	        				}

	    	        				}

	    	        				// BUILD MOTHER VCF //

	    	        				String motherVCF = "./.:.:" + motherQualityDepth;

	    	        				if(motherContents.length >= 15)
	    	        				{
	    	        					if(motherContents.length >= 16)
	    	        						motherAllele = motherContents[15];
	    	        					String consBase = motherContents[0];
		    	        				int reads1 = Integer.parseInt(motherContents[1]);
		    	        				int reads2 = Integer.parseInt(motherContents[2]);
		    	        				String varFreq = motherContents[3];
		    	        				int qual1 = Integer.parseInt(motherContents[6]);
		    	        				int qual2 = Integer.parseInt(motherContents[7]);
		    	        				double pValue = Double.parseDouble(motherContents[8]);
		    	        				int reads1plus = Integer.parseInt(motherContents[11]);
		    	        				int reads1minus = Integer.parseInt(motherContents[12]);
		    	        				int reads2plus = Integer.parseInt(motherContents[13]);
		    	        				int reads2minus = Integer.parseInt(motherContents[14]);

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

		    	        				// mother is wildtype //
		    	        				if(consBase.equals(refBase))
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
		    	        					motherVCF = "0" + "/" + "0";
		    	        					motherVCF += ":" + (int) newLogP + ":" + motherDepth + ":" + motherQualityDepth;
	    	        						motherVCF += ":" + reads1 + ":" + reads2 + ":" + varFreq + ":" + pvalueFormat.format(pValue);
	    	        						motherVCF += ":" + qual1 + ":" + qual2;
	    	        						motherVCF += ":" + reads1plus + ":" + reads1minus + ":" + reads2plus + ":" + reads2minus;

		    	        				}
		    	        				// mother is variant //
		    	        				else if(motherAllele.length() > 0 && !motherAllele.equals("N") && !motherAllele.equals("."))
		    	        				{
		    	        					// Determine how many variant alleles have been seen //

		    	        					int varAlleleNumber = 0;

		    	        					// Determine if we've seen the variant and what its number is ##

			    	        				if(varAlleles.containsKey(motherAllele))
			    	        				{
			    	        					varAlleleNumber = varAlleles.get(motherAllele);
			    	        				}
			    	        				else
			    	        				{
			    	        					// IF no variants yet seen, this is variant allele 1 //
			    	        					varAlleleNumber = varAlleles.size() + 1;
			    	        					varAlleles.put(motherAllele, varAlleleNumber);
			    	        				}

			    	        				if(motherContents.length >= 1)
			    	        				{
			    	        					if(VarScan.isHomozygous(consBase))
			    	        					{
			    	        						motherVCF = varAlleleNumber + "/" + varAlleleNumber;
			    	        					}
			    	        					else
			    	        					{
			    	        						motherVCF = "0" + "/" + varAlleleNumber;
			    	        					}

			    	        					motherVCF += ":" + (int) logP + ":" + motherDepth + ":" + motherQualityDepth;
		    	        						motherVCF += ":" + reads1 + ":" + reads2 + ":" + varFreq + ":" + pvalueFormat.format(pValue);
		    	        						motherVCF += ":" + qual1 + ":" + qual2;
		    	        						motherVCF += ":" + reads1plus + ":" + reads1minus + ":" + reads2plus + ":" + reads2minus;
			    	        				}

		    	        					if(motherAllele.length() > 1)
		    	        						variantType = "Indel";
		    	        				}

	    	        				}

	    	        				// BUILD CHILD VCF //

	    	        				String childVCF = "./.:.:" + childQualityDepth;

	    	        				if(childContents.length >= 15)
	    	        				{
	    	        					if(childContents.length >= 16)
	    	        						childAllele = childContents[15];
	    	        					String consBase = childContents[0];
		    	        				int reads1 = Integer.parseInt(childContents[1]);
		    	        				int reads2 = Integer.parseInt(childContents[2]);
		    	        				String varFreq = childContents[3];
		    	        				int qual1 = Integer.parseInt(childContents[6]);
		    	        				int qual2 = Integer.parseInt(childContents[7]);
		    	        				double pValue = Double.parseDouble(childContents[8]);
		    	        				int reads1plus = Integer.parseInt(childContents[11]);
		    	        				int reads1minus = Integer.parseInt(childContents[12]);
		    	        				int reads2plus = Integer.parseInt(childContents[13]);
		    	        				int reads2minus = Integer.parseInt(childContents[14]);

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

		    	        				// child is wildtype //
		    	        				if(consBase.equals(refBase))
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
		    	        					childVCF = "0" + "/" + "0";
		    	        					childVCF += ":" + (int) newLogP + ":" + childDepth + ":" + childQualityDepth;
	    	        						childVCF += ":" + reads1 + ":" + reads2 + ":" + varFreq + ":" + pvalueFormat.format(pValue);
	    	        						childVCF += ":" + qual1 + ":" + qual2;
	    	        						childVCF += ":" + reads1plus + ":" + reads1minus + ":" + reads2plus + ":" + reads2minus;

		    	        				}
		    	        				// child is variant //
		    	        				else if(childAllele.length() > 0 && !childAllele.equals("N") && !childAllele.equals("."))
		    	        				{
		    	        					// Determine how many variant alleles have been seen //

		    	        					int varAlleleNumber = 0;

		    	        					// Determine if we've seen the variant and what its number is ##

			    	        				if(varAlleles.containsKey(childAllele))
			    	        				{
			    	        					varAlleleNumber = varAlleles.get(childAllele);
			    	        				}
			    	        				else
			    	        				{
			    	        					// IF no variants yet seen, this is variant allele 1 //
			    	        					varAlleleNumber = varAlleles.size() + 1;
			    	        					varAlleles.put(childAllele, varAlleleNumber);
			    	        				}

			    	        				if(childContents.length >= 1)
			    	        				{
			    	        					if(VarScan.isHomozygous(consBase))
			    	        					{
			    	        						childVCF = varAlleleNumber + "/" + varAlleleNumber;
			    	        					}
			    	        					else
			    	        					{
			    	        						childVCF = "0" + "/" + varAlleleNumber;
			    	        					}

			    	        					childVCF += ":" + (int) logP + ":" + childDepth + ":" + childQualityDepth;
		    	        						childVCF += ":" + reads1 + ":" + reads2 + ":" + varFreq + ":" + pvalueFormat.format(pValue);
		    	        						childVCF += ":" + qual1 + ":" + qual2;
		    	        						childVCF += ":" + reads1plus + ":" + reads1minus + ":" + reads2plus + ":" + reads2minus;
			    	        				}

		    	        					if(childAllele.length() > 1)
		    	        						variantType = "Indel";
		    	        				}

	    	        				}


	    	        				// BEGIN BUILDING OUTPUT //


		    	        			// Build an output line //
		    	        			String outLine = refName + "\t" + position + "\t";


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


		    						// Calculate average sample depth //
		    						int avgQualityDepth = (fatherQualityDepth + motherQualityDepth + childQualityDepth) / 3;
		    						String refColumn = "";
		    						String varColumn = "";

		    						// Handle complex positions with multiple alleles including at least one indel //

		    						if(varBases.contains(",") && (varBases.contains("-") || varBases.contains("+")))
		    						{
		    							variantType = "INDEL";
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
		    							variantType = "INDEL";
			    						// INSERTION //
			    						// Ref = ref base; Var = ref base followed by inserted bases //
		    							refColumn = refBase;
			    						varColumn = refBase + varBases.replace("+", "");
			    					}
			    					else if(varBases.startsWith("-"))
			    					{
			    						variantType = "INDEL";
			    						// DELETION //
			    						// Ref = ref base followed by deleted bases; var = ref base //
			    						refColumn = refBase + varBases.replace("-", "");
			    						varColumn = refBase;
			    					}
			    					else
			    					{
			    						// Variant type SNP //
			    						refColumn = refBase;
			    						varColumn = varBases;
			    					}


		    						// Ensure that varColumn does not contain any +/- //
		    						varColumn = varColumn.replace("+", "");
		    						varColumn = varColumn.replace("-", "");

		    						// ADD REF, ALT, FILTER, INFO, and FORMAT FIELDS TO OUTPUT //

		    						outLine += "." + "\t" + refColumn + "\t" + varColumn + "\t.\t";

		    						String filterColumn = "";
		    						if(trioStatus.contains("MIE"))
		    						{
		    							filterColumn = "mendelError";
		    						}
		    						else if (strandFilterStatus.contains("Fail"))
		    						{
		    							filterColumn = "str10";
		    						}
		    						else
		    						{
		    							filterColumn = "PASS";
		    						}

		    						outLine += filterColumn + "\t";
		    						outLine += "ADP=" + avgQualityDepth + ";STATUS="; // + trioStatus;

		    						if(trioStatus.contains("Untransmitted"))
		    						{
		    							outLine += "1";
		    						}
		    						else if(trioStatus.contains("Germline"))
		    						{
		    							outLine += "2";
		    						}
		    						else if(trioStatus.contains("DeNovo"))
		    						{
		    							outLine += "3;DENOVO";
		    						}
		    						else if(trioStatus.contains("MIE"))
		    						{
		    							outLine += "4";
		    						}

		    						outLine += "\t" + "GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR" + "\t";

		    	        			outLine += fatherVCF + "\t" + motherVCF + "\t" + childVCF;
		    	        			//			    	    			outLine += "\t" + fatherCall.replace("\t", ":") + "\t" + motherCall.replace("\t", ":") + "\t" + childCall.replace("\t", ":");

			    	    			// Count this trio status //
			    	    			String statKey = "";
			    	    			if(trioStatus.equals("Reference"))
			    	    			{
			    	    				// No counting or printing these sites //
			    	    			}
			    	    			else
			    	    			{
			    	    				// A variant position... flag it and count the type //
			    	    				numVariantPositions++;
			    	    				variantFlag = true;

			    	    				if(variantType.equals("INDEL"))
			    	    				{
			    	    					numIndelPositions++;
			    	    				}
			    	    				else
			    	    				{
			    	    					numSNPpositions++;
			    	    				}

			    	    				// Also count pass/fail filter statuses //

			    	    				if(strandFilterStatus.contains("Fail"))
			    	    				{
			    	    					numFailStrandFilter++;
			    	    				}
			    	    				else if(trioStatus.equals("MIE"))
			    	    				{
			    	    					numFailMendelFilter++;
			    	    				}
			    	    				else
			    	    				{
			    	    					numVariantsReported++;
			    	    					if(trioStatus.equals("DeNovo"))
			    	    						numVariantsReportedDeNovo++;

				    	    				if(variantType.equals("INDEL"))
				    	    				{
				    	    					numIndelsReported++;
				    	    					if(trioStatus.equals("DeNovo"))
				    	    						numIndelsReportedDeNovo++;
				    	    				}
				    	    				else
				    	    				{
				    	    					numSNPsReported++;
				    	    					if(trioStatus.equals("DeNovo"))
				    	    						numSNPsReportedDeNovo++;
				    	    				}
			    	    				}

			    	    			}

			    	    			// Determine if we should print the output line //
			    	    			if(variantFlag)
			    	    			{
			    	    				if(variantType.equals("SNP"))
			    	    				{
			    	    					outSnp.println(outLine);
			    	    				}
			    	    				else if(variantType.equals("INDEL"))
			    	    				{
			    	    					outIndel.println(outLine);
			    	    				}
			    	    			}
		    	        		}
		    	        		else
		    	        		{
		    	        			// The trioCallContents was less than 4 fields, so that's a problem //
		    	        			System.err.println("No status for " + numBases);
		    	        		}

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
		    					System.err.println("Error: Invalid format for pileup at line " + numBases + "\n" + line + "\n");
		    					return;
	    					}

	    				}
	    			}
	    			catch(Exception e)
	    		    {
	    		    	System.err.println("Parsing Exception on line:\n" + line + "\n" + e.getMessage() + "\n" + e.getLocalizedMessage());
	    		    	e.printStackTrace();
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
				System.err.println(numBasesCovered + " met the coverage requirement of " + minCoverage);
				System.err.println(numVariantPositions + " variant positions (" + numSNPpositions + " SNP, " + numIndelPositions + " indel)");
				System.err.println(numFailStrandFilter + " were failed by the strand-filter");
				System.err.println(numVariantsReported + " variant positions reported (" + numSNPsReported + " SNP, " + numIndelsReported + " indel)");
				System.err.println(numVariantsReportedDeNovo + " de novo mutations reported (" + numSNPsReportedDeNovo + " SNP, " + numIndelsReportedDeNovo + " indel)");

				// Print the status of each trio call //
//				System.err.println(stats.get("Reference") + " called Reference");
				String[] statsKeys = (String[]) stats.keySet().toArray(new String[0]);

				Arrays.sort(statsKeys);

				// Get the total number of reads at this position //
				int totalstats = 0;
				for(String statsKey : statsKeys)
				{
					System.err.println(stats.get(statsKey) + " " + statsKey);
				}

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





	public String callTrio(String refBase, HashMap<String, String> fatherCounts, HashMap<String, String> motherCounts, HashMap<String, String> childCounts, int minReads2, double minVarFreq, int minAvgQual, double pValueThreshold, double minFreqForHom)
	{
		String fatherCall = VarScan.callPosition(refBase, fatherCounts, "CNS", minReads2, minVarFreq, minAvgQual, pValueThreshold, minFreqForHom);
		String motherCall = VarScan.callPosition(refBase, motherCounts, "CNS", minReads2, minVarFreq, minAvgQual, pValueThreshold, minFreqForHom);
		String childCall = VarScan.callPosition(refBase, childCounts, "CNS", minReads2, minVarFreq, minAvgQual, pValueThreshold, minFreqForHom);

		// Determine the father, mother, and child genotypes //
		String trioStatus = "unknown";

		try
		{
			String[] fatherContents = fatherCall.split("\t");
			String[] motherContents = motherCall.split("\t");
			String[] childContents = childCall.split("\t");

			String father = fatherContents[0];
			String mother = motherContents[0];
			String child = childContents[0];
			String fatherAllele = refBase;
			String motherAllele = refBase;
			String childAllele = refBase;

			if(fatherContents.length >= 16)
				fatherAllele = fatherContents[15];

			if(motherContents.length >= 16)
				motherAllele = motherContents[15];

			if(childContents.length >= 16)
				childAllele = childContents[15];

			// Uninteresting case 1: Any Sample called N //
			if(child.equals("N") || father.equals("N") || mother.equals("N"))
			{
				// Missing data, so not sure
				trioStatus = "MissingData";
			}
			// CASE 1: ALL 3 SAMPLES WILDTYPE //
			else if(father.equals(refBase) && mother.equals(refBase) && child.equals(refBase))
			{
				trioStatus = "Reference";
			}
			// CASE 2: DE NOVO MUTATION //
			else if(father.equals(refBase) && mother.equals(refBase) && !child.equals(refBase))
			{
				trioStatus = "DeNovo";

			}
			else
			{
				// CHECK INDIVIDUAL GENOTYPE ALLELES //
				String fatherGt = VarScan.codeToGenotype(father);
				String motherGt = VarScan.codeToGenotype(mother);
				String childGt = VarScan.codeToGenotype(child);

				String father1 = "N";
				String father2 = "N";
				String[] fatherSplit = fatherGt.split("/");
				if(fatherSplit.length == 2)
				{
					father1 = fatherSplit[0];
					father2 = fatherSplit[1];
				}

				String mother1 = "N";
				String mother2 = "N";
				String[] motherSplit = motherGt.split("/");
				if(motherSplit.length == 2)
				{
					mother1 = motherSplit[0];
					mother2 = motherSplit[1];
				}

				String child1 = "N";
				String child2 = "N";
				String[] childSplit = childGt.split("/");
				if(childSplit.length == 2)
				{
					child1 = childSplit[0];
					child2 = childSplit[1];
				}

				if(father1.equals("*"))
					father1 = refBase;
				if(mother1.equals("*"))
					mother1 = refBase;
				if(child1.equals("*"))
					child1 = refBase;

				// CHILD IS VARIANT
				if(!child.equals(refBase))
				{
					if((child1.equals(father1) || child1.equals(father2) || child2.equals(father1) || child2.equals(father2)) && (child1.equals(mother1) || child1.equals(mother2) || child2.equals(mother1) || child2.equals(mother2)))
					{
						trioStatus = "Germline";
					}
					else if(!child.equals(fatherAllele) && !child.equals(motherAllele))
					{
						trioStatus = "MultAlleles";
					}
					else
					{
						trioStatus = "MIE";
					}
				}
				// CHILD IS WILDTYPE
				else if(child.equals(refBase))
				{
					// Can one wildtype allele come from dad and one from mom ? //
					if((father1.equals(refBase) || father2.equals(refBase)) && (mother1.equals(refBase) || mother2.equals(refBase)))
					{
						trioStatus = "Untransmitted";
					}

					else
					{
						trioStatus = "MIE";
					}
				}
			}
			// CASE 3B: GERMLINE VARIANT WITH CHILD HET AND ONE PARENT HOM
			// CASE 3C: CHILD IS WILDTYPE, AND NEITHER PARENT HOMOZYGOUS-VARIANT
			// CASE 3C: CHILD IS WILDTYPE, AND NEITHER PARENT HOMOZYGOUS-VARIANT
			// CASE 4: IMPOSSIBLE HOMOZYGOTE
			// CASE 5: SHOULD BE HET //

			return(fatherCall.replace("\t", ":") + "\t" + motherCall.replace("\t", ":") + "\t" + childCall.replace("\t", ":") + "\t" + trioStatus);

		}
		catch(Exception e)
		{
			System.err.println("Error parsing genotypes: " + e.getMessage() + " local: " + e.getLocalizedMessage());
			e.printStackTrace();

		}

		return("");
	}

}
