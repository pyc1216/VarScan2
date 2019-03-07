/**
 * @(#)Somatic.java
 *
 * Copyright (c) 2009-2010 Daniel C. Koboldt and Washington University in St. Louis
 *
 * COPYRIGHT
 */

package net.sf.varscan;

//Import required packages //

import java.io.*;
import java.text.*;
import java.util.*;
import java.lang.Math.*;

/**
 * A class for determining somatic status of variants from Normal/Tumor pileup files
 *
 * @version	2.3
 *
 * @author Daniel C. Koboldt <dkoboldt@genome.wustl.edu>
 *
 */
public class Somatic {




	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Constructor with two arguments (string[], boolean) expects mpileup input 					  //
	////////////////////////////////////////////////////////////////////////////////////////////////////

	public Somatic(String[] args, boolean isMpileup)
	{
		String usage = "USAGE: java -jar VarScan.jar somatic [normal-tumor.mpileup] [Opt: output] OPTIONS\n" +
			"\tnormal-tumor.pileup - The SAMtools mpileup file for Normal and Tumor BAMs\n" +
			"\toutput - Output base name for SNP and indel output\n" +
			"\nOPTIONS:\n" +
			"\t--output-snp - Output file for SNP calls [output.snp]\n" +
			"\t--output-indel - Output file for indel calls [output.indel]\n" +
			"\t--min-coverage - Minimum coverage in normal and tumor to call variant [8]\n" +
			"\t--min-coverage-normal - Minimum coverage in normal to call somatic [8]\n" +
			"\t--min-coverage-tumor - Minimum coverage in tumor to call somatic [6]\n" +
			"\t--min-avg-qual - Minimum Phred quality to count a base [15]\n" +
			"\t--min-var-freq - Minimum variant frequency to call a heterozygote [0.10]\n" +
			"\t--min-freq-for-hom\tMinimum frequency to call homozygote [0.75]\n" +
			"\t--normal-purity - Estimated purity (non-tumor content) of normal sample [1.00]\n" +
			"\t--tumor-purity - Estimated purity (tumor content) of tumor sample [1.00]\n" +
			"\t--p-value - P-value threshold to call a heterozygote [0.99]\n" +
			"\t--somatic-p-value - P-value threshold to call a somatic site [0.05]\n" +
			"\t--strand-filter - If set to 1, removes variants with >90% strand bias\n" +
			"\t--validation - If set to 1, outputs all compared positions even if non-variant\n" +
			"\t--output-vcf - If set to 1, output VCF instead of VarScan native format\n";

		String vcfHeader = "##fileformat=VCFv4.1";
		vcfHeader += "\n" + "##source=VarScan2";
		vcfHeader += "\n" + "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total depth of quality bases\">";
		vcfHeader += "\n" + "##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description=\"Indicates if record is a somatic mutation\">";
		vcfHeader += "\n" + "##INFO=<ID=SS,Number=1,Type=String,Description=\"Somatic status of variant (0=Reference,1=Germline,2=Somatic,3=LOH, or 5=Unknown)\">";
		vcfHeader += "\n" + "##INFO=<ID=SSC,Number=1,Type=String,Description=\"Somatic score in Phred scale (0-255) derived from somatic p-value\">";
		vcfHeader += "\n" + "##INFO=<ID=GPV,Number=1,Type=Float,Description=\"Fisher's Exact Test P-value of tumor+normal versus no variant for Germline calls\">";
		vcfHeader += "\n" + "##INFO=<ID=SPV,Number=1,Type=Float,Description=\"Fisher's Exact Test P-value of tumor versus normal for Somatic/LOH calls\">";
		vcfHeader += "\n" + "##FILTER=<ID=str10,Description=\"Less than 10% or more than 90% of variant supporting reads on one strand\">";
		vcfHeader += "\n" + "##FILTER=<ID=indelError,Description=\"Likely artifact due to indel reads at this position\">";
		vcfHeader += "\n" + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">";
		vcfHeader += "\n" + "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">";
		vcfHeader += "\n" + "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">";
		vcfHeader += "\n" + "##FORMAT=<ID=RD,Number=1,Type=Integer,Description=\"Depth of reference-supporting bases (reads1)\">";
		vcfHeader += "\n" + "##FORMAT=<ID=AD,Number=1,Type=Integer,Description=\"Depth of variant-supporting bases (reads2)\">";
		vcfHeader += "\n" + "##FORMAT=<ID=FREQ,Number=1,Type=String,Description=\"Variant allele frequency\">";
		vcfHeader += "\n" + "##FORMAT=<ID=DP4,Number=1,Type=String,Description=\"Strand read counts: ref/fwd, ref/rev, var/fwd, var/rev\">";
		vcfHeader += "\n" + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\tTUMOR";

		// Set parameter defaults //

		HashMap<String, String> params = VarScan.getParams(args);

		// Set up formatting for p-values //
		DecimalFormat pvalueFormat = new DecimalFormat("0.####E0");
	    DecimalFormat oneDigit = new DecimalFormat("#0.0");
	    DecimalFormat threeDigits = new DecimalFormat("#0.000");

		// Establish output file names //
		String outputName = "output";
		String outputSnp = "";
		String outputIndel = "";
		String outputCopy = "";

		if(args.length >= 3 && !args[2].startsWith("-"))
		{
			outputName = args[2];
			outputSnp = outputName + ".snp";
			outputIndel = outputName + ".indel";
		}

//		Set parameter defaults //

		int minCoverage = 	8;
		int minCoverageNormal = 8;
		int minCoverageTumor = 6;
		int minReads2 = 	2;
		int minStrands2 = 	1;
		int minAvgQual = 	15;
		double normalPurity = 1.00;
		double tumorPurity = 1.00;
		double dataRatio = 1.00;
		double minVarFreq = 0.20;
		double pValueThreshold = 0.99;
		double somaticPvalue = 0.05; //1.0e-04;
		double minFreqForHom = 0.75;
		boolean doStrandFilter = true;

		// Try adjusting any provided parameters based on user inut //
		try
		{
			if(params.containsKey("output-snp"))
				 outputSnp = params.get("output-snp");

			if(params.containsKey("output-indel"))
				 outputIndel = params.get("output-indel");

			if(params.containsKey("min-coverage"))
			{
				 minCoverage = Integer.parseInt(params.get("min-coverage"));
				 minCoverageNormal = minCoverage;
				 minCoverageTumor = minCoverage;
			}

			if(params.containsKey("min-coverage-normal"))
				 minCoverageNormal = Integer.parseInt(params.get("min-coverage-normal"));

			if(params.containsKey("min-coverage-tumor"))
				 minCoverageTumor = Integer.parseInt(params.get("min-coverage-tumor"));

			if(params.containsKey("min-reads2"))
				 minReads2 = Integer.parseInt(params.get("min-reads2"));

			if(params.containsKey("min-strands2"))
				 minStrands2 = Integer.parseInt(params.get("min-strands2"));

			if(params.containsKey("min-var-freq"))
				 minVarFreq = Double.parseDouble(params.get("min-var-freq"));

			if(params.containsKey("min-freq-for-hom"))
				 minFreqForHom = Double.parseDouble(params.get("min-freq-for-hom"));

			if(params.containsKey("min-avg-qual"))
				 minAvgQual = Integer.parseInt(params.get("min-avg-qual"));

			if(params.containsKey("p-value"))
				 pValueThreshold = Double.parseDouble(params.get("p-value"));

			if(params.containsKey("somatic-p-value"))
				 somaticPvalue = Double.parseDouble(params.get("somatic-p-value"));

			if(params.containsKey("data-ratio"))
				 dataRatio = Double.parseDouble(params.get("data-ratio"));

			if(params.containsKey("normal-purity"))
			{
				 normalPurity = Double.parseDouble(params.get("normal-purity"));
				 if(normalPurity > 1)
					 normalPurity = normalPurity / 100.00;
			}

			if(params.containsKey("tumor-purity"))
			{
				 tumorPurity = Double.parseDouble(params.get("tumor-purity"));
				 if(tumorPurity > 1)
					 tumorPurity = normalPurity / 100.00;
			}

			if(params.containsKey("strand-filter"))
			{
				int filter = Integer.parseInt(params.get("strand-filter"));
				if(filter > 0)
					doStrandFilter = true;
				else
					doStrandFilter = false;
			}

//			System.err.println("Min coverage:\t" + minCoverage);
			System.err.println("Min coverage:\t" + minCoverageNormal + "x for Normal, " + minCoverageTumor + "x for Tumor");
			System.err.println("Min reads2:\t" + minReads2);
			System.err.println("Min strands2:\t" + minStrands2);
			System.err.println("Min var freq:\t" + minVarFreq);
			System.err.println("Min freq for hom:\t" + minFreqForHom);
			System.err.println("Normal purity:\t" + normalPurity);
			System.err.println("Tumor purity:\t" + tumorPurity);
			System.err.println("Min avg qual:\t" + minAvgQual);
			System.err.println("P-value thresh:\t" + pValueThreshold);
			System.err.println("Somatic p-value:\t" + somaticPvalue);
			if(params.containsKey("validation"))
				System.err.println("Validation mode: on");

		}
		catch(Exception e)
		{
	    	System.err.println("Input Parameter Threw Exception: " + e.getLocalizedMessage());
	    	e.printStackTrace(System.err);
	    	System.exit(1);
		}

		// Print usage if -h or --help invoked //
		if(params.containsKey("help") || params.containsKey("h"))
		{
			System.err.println(usage);
			return;
		}

		// Check for correct input //

		if(outputSnp.length() == 0 || outputIndel.length() == 0)
		{
			System.err.println("Please provide an output basename or SNP/indel output files!");
			 System.err.println(usage);
			 System.exit(1);
		}

		// Statistics counters //
		long tumorPositions = 0;
		long sharedPositions = 0;
		long comparedPositions = 0;
		long calledReference = 0;
		long indelFilter = 0;
		long strandFilter = 0;
		long calledGermline = 0;
		long calledLOH = 0;
		long calledSomatic = 0;
		long calledUnknown = 0;
		long calledVariant = 0;


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

	    	if(in != null && in.ready())
	    	{
	    		// Declare some file-parsing variables //
			    String lineNormal;
			    String lineTumor;
			    String chromNormal = "";
			    String chromTumor = "";
			    String refBase = "";
			    int posNormal = 0;
			    int posTumor = 0;

				// Declare output file //
		 	 	PrintStream outSnp = null; // declare a print stream object for SNPs
		 	 	PrintStream outIndel = null; // declare a print stream object for Indels
		 	 	PrintStream outValidation = null; // declare a print stream object for both for validation
		 	 	PrintStream outCopyNumber = null; // declare a print stream object for both for validation

		 	 	if(params.containsKey("output-vcf"))
		 	 	{
		 	 		if(!outputSnp.contains(".vcf"))
		 	 			outputSnp += ".vcf";
		 	 		if(!outputIndel.contains(".vcf"))
		 	 			outputIndel += ".vcf";
		 	 	}
		 		outSnp = new PrintStream( new FileOutputStream(outputSnp) );
		 		outIndel = new PrintStream( new FileOutputStream(outputIndel) );

		 		if(!params.containsKey("no-headers") && !params.containsKey("output-vcf"))
		 		{
		 			outSnp.println("chrom\tposition\tref\tvar\tnormal_reads1\tnormal_reads2\tnormal_var_freq\tnormal_gt\ttumor_reads1\ttumor_reads2\ttumor_var_freq\ttumor_gt\tsomatic_status\tvariant_p_value\tsomatic_p_value\ttumor_reads1_plus\ttumor_reads1_minus\ttumor_reads2_plus\ttumor_reads2_minus\tnormal_reads1_plus\tnormal_reads1_minus\tnormal_reads2_plus\tnormal_reads2_minus");
		 			outIndel.println("chrom\tposition\tref\tvar\tnormal_reads1\tnormal_reads2\tnormal_var_freq\tnormal_gt\ttumor_reads1\ttumor_reads2\ttumor_var_freq\ttumor_gt\tsomatic_status\tvariant_p_value\tsomatic_p_value\ttumor_reads1_plus\ttumor_reads1_minus\ttumor_reads2_plus\ttumor_reads2_minus\tnormal_reads1_plus\tnormal_reads1_minus\tnormal_reads2_plus\tnormal_reads2_minus");
		 		}

		 		if(params.containsKey("output-vcf"))
				{
					// Output VCF Header //
					outSnp.println(vcfHeader);
					outIndel.println(vcfHeader);
				}

		 		if(params.containsKey("validation"))
		 		{
			 		outValidation = new PrintStream( new FileOutputStream(outputName + ".validation") );
			 		if(!params.containsKey("no-headers") && !params.containsKey("output-vcf"))
			 			outValidation.println("chrom\tposition\tref\tvar\tnormal_reads1\tnormal_reads2\tnormal_var_freq\tnormal_gt\ttumor_reads1\ttumor_reads2\ttumor_var_freq\ttumor_gt\tsomatic_status\tvariant_p_value\tsomatic_p_value\ttumor_reads1_plus\ttumor_reads1_minus\ttumor_reads2_plus\ttumor_reads2_minus");
			 		if(params.containsKey("output-vcf"))
					{
						// Output VCF Header //
						outValidation.println(vcfHeader);
					}
		 		}

	    		// Parse the infile line by line //
	    		System.err.println("Reading mpileup input...");
	    		int numParsingExceptions = 0;

	    		while ((line = in.readLine()) != null)
	    		{

	    			// Begin try-catch for line parsing //

	    			try
	    			{
	    				String[] lineContents = line.split("\t", -1);

	    				// Verify expected pileup format //

	    				if(lineContents.length > 5 && lineContents[0].length() > 0 && lineContents[1].length() > 0 && lineContents[2].length() > 0 && lineContents[3].length() > 0)
	    				{
	    					sharedPositions++;

	    					// Parse common fields from line //
	    					String refName = lineContents[0];
	    					int position = Integer.parseInt(lineContents[1]);
	    					refBase = lineContents[2].toUpperCase();

	    					chromNormal = refName;
	    					chromTumor = refName;
	    					posNormal = position;
	    					posTumor = position;

	    					// Parse normal, which should be first sample //
	    					int normalOffset = 3;
	    					int pileupDepthNormal = Integer.parseInt(lineContents[normalOffset]);
		    	        	String normalBases = lineContents[normalOffset + 1];
		    	        	String normalQualities = lineContents[normalOffset + 2];

	    					// Parse tumor, which should be second sample //
	    					int tumorOffset = 6;
	    					int pileupDepthTumor = Integer.parseInt(lineContents[tumorOffset]);
		    	        	String tumorBases = lineContents[tumorOffset + 1];
		    	        	String tumorQualities = lineContents[tumorOffset + 2];

	    					lineNormal = refName + "\t" + position + "\t" + refBase + "\t" + pileupDepthNormal + "\t" + normalBases + "\t" + normalQualities;
	    					lineTumor = refName + "\t" + position + "\t" + refBase + "\t" + pileupDepthTumor + "\t" + tumorBases + "\t" + tumorQualities;

	    					String compareResult = comparePositions(lineNormal, lineTumor, minCoverage, minReads2, minVarFreq, minAvgQual, pValueThreshold, somaticPvalue, minFreqForHom, normalPurity, tumorPurity);

			    			if(compareResult.length() > 0)
			    			{
			    				// Get the alleles to determine type //
			    				String[] compareContents = compareResult.split("\t");
				    			String allele1 = compareContents[0];
				    			String allele2 = compareContents[1];

				    			double strandedness1 = 0.50;
				    			double strandedness2 = 0.50;
				    			double strandednessDiff = 0.00;

				    			if(compareContents.length >= 17)
				    			{
				    				try
				    				{
				    					int tumorReads1plus = Integer.parseInt(compareContents[13]);
				    					int tumorReads1minus = Integer.parseInt(compareContents[14]);
				    					int tumorReads2plus = Integer.parseInt(compareContents[15]);
				    					int tumorReads2minus = Integer.parseInt(compareContents[16]);

				    					if(tumorReads1plus > 0 || tumorReads1minus > 0)
				    					{
				    						strandedness1 = (double) tumorReads1plus / (double) (tumorReads1plus + tumorReads1minus);
				    					}

				    					if(tumorReads2plus > 0 || tumorReads2minus > 0)
				    					{
				    						strandedness2 = (double) tumorReads2plus / (double) (tumorReads2plus + tumorReads2minus);
				    						if(tumorReads1plus > 0 || tumorReads1minus > 0)
				    						{
				    							strandednessDiff = java.lang.Math.abs(strandedness1 - strandedness2);
				    						}
				    					}
				    				}
				    				catch(Exception e)
				    				{
				    					// Exception parsing info from compareResult //
				    				}
				    			}

			    				//stats.put("comparedPositions", (stats.get("comparedPositions") + 1));
				    			comparedPositions++;

			    				if(params.containsKey("verbose") && !compareResult.contains("Reference"))
			    					System.err.println(chromNormal + "\t" + posNormal + "\t" + compareResult);

			    				// If VCF format specified, supply it //

			    				if(params.containsKey("output-vcf"))
			    				{
			    					int normalReads1 = Integer.parseInt(compareContents[2]);
			    					int normalReads2 = Integer.parseInt(compareContents[3]);
			    					String normalFreq = compareContents[4];
			    					String normalCall = compareContents[5];
			    					int tumorReads1 = Integer.parseInt(compareContents[6]);
			    					int tumorReads2 = Integer.parseInt(compareContents[7]);
			    					String tumorFreq = compareContents[8];
			    					String tumorCall = compareContents[9];
			    					String somStatus = compareContents[10];
			    					Double germlineP = Double.parseDouble(compareContents[11]);
			    					Double somaticP = Double.parseDouble(compareContents[12]);

			    					int totalDepth = pileupDepthNormal + pileupDepthTumor;

			    					if(allele2.startsWith("+"))
			    					{
			    						// INSERTION //
			    						// Ref = ref base; Var = ref base followed by inserted bases //
			    						String varColumn = allele1 + allele2.replace("+", "");
			    						compareResult = "." + "\t" + allele1 + "\t" + varColumn + "\t" + ".";
			    					}
			    					else if(allele2.startsWith("-"))
			    					{
			    						// DELETION //
			    						// Ref = ref base followed by deleted bases; var = ref base //
			    						String refColumn = allele1 + allele2.replace("-", "");
			    						compareResult = "." + "\t" + refColumn + "\t" + allele1 + "\t" + ".";
			    					}
			    					else
			    					{
				    					compareResult = "." + "\t" + allele1 + "\t" + allele2 + "\t" + ".";
			    					}


			    					// Decide on filter field //
			    					if(doStrandFilter && strandednessDiff > 0.10 && (strandedness2 < 0.10 || strandedness2 > 0.90))
			    					{
			    						compareResult += "\t" + "str10";
			    					}
			    					else if(somStatus.equals("IndelFilter"))
			    					{
			    						compareResult += "\t" + "indelError";
			    					}
			    					else
			    					{
			    						compareResult += "\t" + "PASS";
			    					}

			    					// Determine somatic status id and score //
			    					int ssCode = 0;
			    					double somScore = 0;

			    					if(somStatus.equals("Reference"))
			    					{
			    						// Wildtype //
			    						ssCode = 0;
			    						calledReference++;
			    					}
			    					else if(somStatus.equals("Germline"))
			    					{
			    						// Germline //
			    						ssCode = 1;
			    						calledGermline++;
			    						if(somaticP == 0)
			    						{
			    							somScore = 0;
			    						}
			    						else
			    						{
			    							somScore = 0 - (10 * java.lang.Math.log10(somaticP));
			    						}
			    					}
			    					else if(somStatus.equals("Somatic"))
			    					{
			    						// Somatic //
			    						ssCode = 2;
			    						calledSomatic++;
			    						if(somaticP == 0)
			    						{
			    							somScore = 255;
			    						}
			    						else
			    						{
			    							somScore = 0 - (10 * java.lang.Math.log10(somaticP));
			    						}
			    					}
			    					else if(somStatus.equals("LOH"))
			    					{
			    						// LOH //
			    						ssCode = 3;
			    						calledLOH++;
			    						if(somaticP == 0)
			    						{
			    							somScore = 255;
			    						}
			    						else
			    						{
			    							somScore = 0 - (10 * java.lang.Math.log10(somaticP));
			    						}
			    					}
			    					else
			    					{
			    						// Unknown //
			    						calledUnknown++;
			    						ssCode = 5;
			    					}

			    					// Adjust somatic score //
			    					if(somScore > 255)
			    						somScore = 255;

			    					// Print the info field //

			    					compareResult += "\t" + "DP=" + totalDepth;
			    					if(somStatus.equals("Somatic"))
			    						compareResult += ";SOMATIC";
			    					compareResult += ";" + "SS=" + ssCode;
			    					compareResult += ";" + "SSC=" + (int) somScore;
			    					compareResult += ";" + "GPV=" + pvalueFormat.format(germlineP);
			    					compareResult += ";" + "SPV=" + pvalueFormat.format(somaticP);

			    					// Print the format field //

			    					String tumorDP4 = "";
			    					String normalDP4 = "";

					    			if(compareContents.length >= 17)
					    			{
					    				try
					    				{
					    					tumorDP4 = compareContents[13] + "," + compareContents[14] + "," + compareContents[15] + "," + compareContents[16];
					    					normalDP4 = compareContents[17] + "," + compareContents[18] + "," + compareContents[19] + "," + compareContents[20];
					    				}
					    				catch(Exception e)
					    				{
					    					// Exception parsing info from compareResult //
					    					tumorDP4 = "";
					    					normalDP4 = "";
					    				}
					    			}

					    			if(tumorDP4.length() > 0)
					    				compareResult += "\tGT:GQ:DP:RD:AD:FREQ:DP4";
					    			else
					    				compareResult += "\tGT:GQ:DP:RD:AD:FREQ";

			    					// Determine normal genotype //
			    					String normalGt = ".";
			    					String tumorGt = ".";
			    					if(normalCall.equals(refBase))
			    					{
			    						normalGt = "0/0";
			    					}
			    					else if(VarScan.isHeterozygous(normalCall))
			    					{
			    						normalGt = "0/1";
			    					}
			    					else
			    					{
			    						normalGt = "1/1";
			    					}

			    					if(tumorCall.equals(refBase))
			    					{
			    						tumorGt = "0/0";
			    					}
			    					else if(VarScan.isHeterozygous(tumorCall))
			    					{
			    						tumorGt = "0/1";
			    					}
			    					else
			    					{
			    						tumorGt = "1/1";
			    					}

			    					if(tumorDP4.length() > 0)
			    					{
				    					compareResult += "\t" + normalGt + ":.:" + pileupDepthNormal + ":" + normalReads1 + ":" + normalReads2 + ":" + normalFreq + ":" + normalDP4;
				    					compareResult += "\t" + tumorGt + ":.:" + pileupDepthTumor + ":" + tumorReads1 + ":" + tumorReads2 + ":" + tumorFreq + ":" + tumorDP4;
			    					}
			    					else
			    					{
				    					compareResult += "\t" + normalGt + ":.:" + pileupDepthNormal + ":" + normalReads1 + ":" + normalReads2 + ":" + normalFreq;
				    					compareResult += "\t" + tumorGt + ":.:" + pileupDepthTumor + ":" + tumorReads1 + ":" + tumorReads2 + ":" + tumorFreq;
			    					}
			    				}

			    				// Print to master file for validation //

			    				if(params.containsKey("validation"))
			    				{
			    					outValidation.println(chromNormal + "\t" + posNormal + "\t" + compareResult);
			    				}

			    				if(!params.containsKey("validation") && (compareResult.contains("Reference") || compareResult.contains("SS=0")  || compareResult.contains("Filter")))
			    				{
			    					// Don't print reference/indelfilter positions unless doing validation //
			    				}
			    				else if(doStrandFilter && strandednessDiff > 0.10 && (strandedness2 < 0.10 || strandedness2 > 0.90))
			    				{
			    					// If filter is set, ignore variants that are supported largely by one strand //
			    					if(!params.containsKey("output-vcf"))
			    						compareResult = "StrandFilter";
			    				}
			    				else if(allele1.contains("-") || allele1.contains("+") || allele2.contains("-") || allele2.contains("+"))//if(compareResult.contains("INS") || compareResult.contains("DEL"))
			    				{
			    					outIndel.println(chromNormal + "\t" + posNormal + "\t" + compareResult);
			    				}
			    				else
			    				{
				    				outSnp.println(chromNormal + "\t" + posNormal + "\t" + compareResult);
			    				}
			    			}
			    			else
			    			{
//			    				System.err.println("Uncalled" + chromNormal + "\t" + posNormal + "\t" + compareResult);
			    			}

			    			if(compareResult.contains("Reference"))
			    				calledReference++; //stats.put("calledReference", (stats.get("calledReference") + 1));
			    			else if(compareResult.contains("IndelFilter"))
			    				indelFilter++;	//stats.put("indelFilter", (stats.get("indelFilter") + 1));
			    			else if(compareResult.contains("StrandFilter"))
			    				strandFilter++;
			    			else if(compareResult.contains("Germline"))
			    				calledGermline++;	//stats.put("calledGermline", (stats.get("calledGermline") + 1));
			    			else if(compareResult.contains("Somatic"))
			    				calledSomatic++;	//stats.put("calledSomatic", (stats.get("calledSomatic") + 1));
			    			else if(compareResult.contains("LOH"))
			    				calledLOH++;	//stats.put("calledLOH", (stats.get("calledLOH") + 1));
			    			else if(compareResult.contains("Unknown"))
			    				calledUnknown++;	//stats.put("calledUnknown", (stats.get("calledUnknown") + 1));
			    			else if(compareResult.contains("Variant"))
			    				calledVariant++;	//stats.put("calledVariant", (stats.get("calledVariant") + 1));



	    				}
	    				else
	    				{
	    					System.err.println("Error: Invalid format or not enough samples in mpileup: " + line + "\n");
	    					return;
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

	    		// Close input/output files //
	    		in.close();
			    outSnp.close();
			    outIndel.close();

			    System.err.println(sharedPositions + " positions in mpileup file"); //stats.get("sharedPositions")
			    System.err.println(comparedPositions + " had sufficient coverage for comparison"); //stats.get("comparedPositions")
			    System.err.println(calledReference + " were called Reference"); //stats.get("calledReference")
			    System.err.println(indelFilter + " were mixed SNP-indel calls and filtered");
			    if(doStrandFilter)
			    	System.err.println(strandFilter + " were removed by the strand filter");
			    System.err.println(calledGermline + " were called Germline");
			    System.err.println(calledLOH + " were called LOH");
			    System.err.println(calledSomatic + " were called Somatic");
			    System.err.println(calledUnknown + " were called Unknown");
			    System.err.println(calledVariant + " were called Variant");
	    	}
	    	else
	    	{
	    		System.err.println("Input file never ready for parsing (maybe due to file I/O)...");
	    		System.exit(10);
	    	}
		}
		catch (IOException e)
		{
			System.err.println("File Parsing Exception: " + e.getLocalizedMessage());
			e.printStackTrace(System.err);
			System.exit(11);
		}

	}

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Constructor with one argument (string[]) expects two-pileup input 					  //
	////////////////////////////////////////////////////////////////////////////////////////////////////

	public Somatic(String[] args)
	{
		String usage = "USAGE: VarScan somatic [normal_pileup] [tumor_pileup] [Opt: output] OPTIONS\n" +
			"\tnormal_pileup - The SAMtools pileup file for Normal\n" +
			"\ttumor_pileup - The SAMtools pileup file for Tumor\n" +
			"\toutput - Output base name for SNP and indel output\n" +
			"\nOPTIONS:\n" +
			"\t--output-snp - Output file for SNP calls [output.snp]\n" +
			"\t--output-indel - Output file for indel calls [output.indel]\n" +
			"\t--min-coverage - Minimum coverage in normal and tumor to call variant [8]\n" +
			"\t--min-coverage-normal - Minimum coverage in normal to call somatic [8]\n" +
			"\t--min-coverage-tumor - Minimum coverage in tumor to call somatic [6]\n" +
			"\t--min-var-freq - Minimum variant frequency to call a heterozygote [0.10]\n" +
			"\t--min-freq-for-hom\tMinimum frequency to call homozygote [0.75]\n" +
			"\t--normal-purity - Estimated purity (non-tumor content) of normal sample [1.00]\n" +
			"\t--tumor-purity - Estimated purity (tumor content) of tumor sample [1.00]\n" +
			"\t--p-value - P-value threshold to call a heterozygote [0.99]\n" +
			"\t--somatic-p-value - P-value threshold to call a somatic site [0.05]\n" +
			"\t--strand-filter - If set to 1, removes variants with >90% strand bias [0]\n" +
			"\t--validation - If set to 1, outputs all compared positions even if non-variant\n" +
			"\t--output-vcf - If set to 1, output VCF instead of VarScan native format\n";

		String vcfHeader = "##fileformat=VCFv4.1";
		vcfHeader += "\n" + "##source=VarScan2";
		vcfHeader += "\n" + "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total depth of quality bases\">";
		vcfHeader += "\n" + "##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description=\"Indicates if record is a somatic mutation\">";
		vcfHeader += "\n" + "##INFO=<ID=SS,Number=1,Type=String,Description=\"Somatic status of variant (0=Reference,1=Germline,2=Somatic,3=LOH, or 5=Unknown)\">";
		vcfHeader += "\n" + "##INFO=<ID=SSC,Number=1,Type=String,Description=\"Somatic score in Phred scale (0-255) derived from somatic p-value\">";
		vcfHeader += "\n" + "##INFO=<ID=GPV,Number=1,Type=Float,Description=\"Fisher's Exact Test P-value of tumor+normal versus no variant for Germline calls\">";
		vcfHeader += "\n" + "##INFO=<ID=SPV,Number=1,Type=Float,Description=\"Fisher's Exact Test P-value of tumor versus normal for Somatic/LOH calls\">";
		vcfHeader += "\n" + "##FILTER=<ID=str10,Description=\"Less than 10% or more than 90% of variant supporting reads on one strand\">";
		vcfHeader += "\n" + "##FILTER=<ID=indelError,Description=\"Likely artifact due to indel reads at this position\">";
		vcfHeader += "\n" + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">";
		vcfHeader += "\n" + "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">";
		vcfHeader += "\n" + "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">";
		vcfHeader += "\n" + "##FORMAT=<ID=RD,Number=1,Type=Integer,Description=\"Depth of reference-supporting bases (reads1)\">";
		vcfHeader += "\n" + "##FORMAT=<ID=AD,Number=1,Type=Integer,Description=\"Depth of variant-supporting bases (reads2)\">";
		vcfHeader += "\n" + "##FORMAT=<ID=FREQ,Number=1,Type=String,Description=\"Variant allele frequency\">";
		vcfHeader += "\n" + "##FORMAT=<ID=DP4,Number=1,Type=String,Description=\"Strand read counts: ref/fwd, ref/rev, var/fwd, var/rev\">";
		vcfHeader += "\n" + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\tTUMOR";

		if(args.length < 3)
		{
			 System.err.println(usage);
			 System.exit(1);
		}

		// Get the required arguments //
		String normalPileupFile = args[1];
		String tumorPileupFile = args[2];

		String outputName = "output";
		String outputSnp = "";
		String outputIndel = "";
		String outputCopy = "";

		if(args.length >= 4 && !args[3].startsWith("-"))
		{
			outputName = args[3];
			outputSnp = outputName + ".snp";
			outputIndel = outputName + ".indel";
		}


		System.err.println("Normal Pileup: " + normalPileupFile);
		System.err.println("Tumor Pileup: " + tumorPileupFile);

		//	Set parameter defaults //

		int minCoverage = 	8;
		int minCoverageNormal = 8;
		int minCoverageTumor = 6;
		int minReads2 = 	2;
		int minStrands2 = 	1;
		int minAvgQual = 	15;
		double normalPurity = 1.00;
		double tumorPurity = 1.00;
		double dataRatio = 1.00;
		double minVarFreq = 0.20;
		double pValueThreshold = 0.99;
		double somaticPvalue = 0.05; //1.0e-04;
		double minFreqForHom = 0.75;

		// Parse command-line parameters //
		HashMap<String, String> params = VarScan.getParams(args);

		// Try adjusting any provided parameters based on user inut //
		try
		{
			if(params.containsKey("output-snp"))
				 outputSnp = params.get("output-snp");

			if(params.containsKey("output-indel"))
				 outputIndel = params.get("output-indel");

			if(params.containsKey("min-coverage"))
			{
				 minCoverage = Integer.parseInt(params.get("min-coverage"));
				 minCoverageNormal = minCoverage;
				 minCoverageTumor = minCoverage;
			}

			if(params.containsKey("min-coverage-normal"))
				 minCoverageNormal = Integer.parseInt(params.get("min-coverage-normal"));

			if(params.containsKey("min-coverage-tumor"))
				 minCoverageTumor = Integer.parseInt(params.get("min-coverage-tumor"));

			if(params.containsKey("min-reads2"))
				 minReads2 = Integer.parseInt(params.get("min-reads2"));

			if(params.containsKey("min-strands2"))
				 minStrands2 = Integer.parseInt(params.get("min-strands2"));

			if(params.containsKey("min-var-freq"))
				 minVarFreq = Double.parseDouble(params.get("min-var-freq"));

			if(params.containsKey("min-freq-for-hom"))
				 minFreqForHom = Double.parseDouble(params.get("min-freq-for-hom"));

			if(params.containsKey("min-avg-qual"))
				 minAvgQual = Integer.parseInt(params.get("min-avg-qual"));

			if(params.containsKey("p-value"))
				 pValueThreshold = Double.parseDouble(params.get("p-value"));

			if(params.containsKey("somatic-p-value"))
				 somaticPvalue = Double.parseDouble(params.get("somatic-p-value"));

			if(params.containsKey("data-ratio"))
				 dataRatio = Double.parseDouble(params.get("data-ratio"));

			if(params.containsKey("normal-purity"))
			{
				 normalPurity = Double.parseDouble(params.get("normal-purity"));
				 if(normalPurity > 1)
					 normalPurity = normalPurity / 100.00;
			}

			if(params.containsKey("tumor-purity"))
			{
				 tumorPurity = Double.parseDouble(params.get("tumor-purity"));
				 if(tumorPurity > 1)
					 tumorPurity = normalPurity / 100.00;
			}

//			System.err.println("Min coverage:\t" + minCoverage);
			System.err.println("Min coverage:\t" + minCoverageNormal + "x for Normal, " + minCoverageTumor + "x for Tumor");
			System.err.println("Min reads2:\t" + minReads2);
			System.err.println("Min strands2:\t" + minStrands2);
			System.err.println("Min var freq:\t" + minVarFreq);
			System.err.println("Min freq for hom:\t" + minFreqForHom);
			System.err.println("Normal purity:\t" + normalPurity);
			System.err.println("Tumor purity:\t" + tumorPurity);
			System.err.println("Min avg qual:\t" + minAvgQual);
			System.err.println("P-value thresh:\t" + pValueThreshold);
			System.err.println("Somatic p-value:\t" + somaticPvalue);
			if(params.containsKey("validation"))
				System.err.println("Validation mode: on");

		}
		catch(Exception e)
		{
	    	System.err.println("Input Parameter Threw Exception: " + e.getLocalizedMessage());
	    	e.printStackTrace(System.err);
	    	System.exit(1);
		}

		// Print usage if -h or --help invoked //
		if(params.containsKey("help") || params.containsKey("h"))
		{
			System.err.println(usage);
			return;
		}

		// Check for correct input //

		if(args.length < 3 && (outputSnp.length() == 0 || outputIndel.length() == 0))
		{
			System.err.println("Please provide SNP and Indel output files!");
			 System.err.println(usage);
			 System.exit(1);
		}


		// Statistics counters //
		long tumorPositions = 0;
		long sharedPositions = 0;
		long comparedPositions = 0;
		long calledReference = 0;
		long indelFilter = 0;
		long strandFilter = 0;
		long calledGermline = 0;
		long calledLOH = 0;
		long calledSomatic = 0;
		long calledUnknown = 0;
		long calledVariant = 0;
		DecimalFormat pvalueFormat = new DecimalFormat("0.####E0");

		try
		{
			// Declare output file //
	 	 	PrintStream outSnp = null; // declare a print stream object for SNPs
	 	 	PrintStream outIndel = null; // declare a print stream object for Indels
	 	 	PrintStream outValidation = null; // declare a print stream object for both for validation
	 	 	PrintStream outCopyNumber = null; // declare a print stream object for both for validation

	 	 	if(params.containsKey("output-vcf"))
	 	 	{
	 	 		if(!outputSnp.contains(".vcf"))
	 	 			outputSnp += ".vcf";
	 	 		if(!outputIndel.contains(".vcf"))
	 	 			outputIndel += ".vcf";
	 	 	}

	 		outSnp = new PrintStream( new FileOutputStream(outputSnp) );
	 		outIndel = new PrintStream( new FileOutputStream(outputIndel) );
	 		if(!params.containsKey("no-headers") && !params.containsKey("output-vcf"))
	 		{
	 			outSnp.println("chrom\tposition\tref\tvar\tnormal_reads1\tnormal_reads2\tnormal_var_freq\tnormal_gt\ttumor_reads1\ttumor_reads2\ttumor_var_freq\ttumor_gt\tsomatic_status\tvariant_p_value\tsomatic_p_value\ttumor_reads1_plus\ttumor_reads1_minus\ttumor_reads2_plus\ttumor_reads2_minus\tnormal_reads1_plus\tnormal_reads1_minus\tnormal_reads2_plus\tnormal_reads2_minus");
	 			outIndel.println("chrom\tposition\tref\tvar\tnormal_reads1\tnormal_reads2\tnormal_var_freq\tnormal_gt\ttumor_reads1\ttumor_reads2\ttumor_var_freq\ttumor_gt\tsomatic_status\tvariant_p_value\tsomatic_p_value\ttumor_reads1_plus\ttumor_reads1_minus\ttumor_reads2_plus\ttumor_reads2_minus\tnormal_reads1_plus\tnormal_reads1_minus\tnormal_reads2_plus\tnormal_reads2_minus");
	 		}
			if(params.containsKey("output-vcf"))
			{
				// Output VCF Header //
				outSnp.println(vcfHeader);
				outIndel.println(vcfHeader);
			}

	 		if(params.containsKey("validation"))
	 		{
		 		outValidation = new PrintStream( new FileOutputStream(outputName + ".validation") );
		 		if(!params.containsKey("no-headers") && !params.containsKey("output-vcf"))
		 			outValidation.println("chrom\tposition\tref\tvar\tnormal_reads1\tnormal_reads2\tnormal_var_freq\tnormal_gt\ttumor_reads1\ttumor_reads2\ttumor_var_freq\ttumor_gt\tsomatic_status\tvariant_p_value\tsomatic_p_value\ttumor_reads1_plus\ttumor_reads1_minus\ttumor_reads2_plus\ttumor_reads2_minus\tnormal_reads1_plus\tnormal_reads1_minus\tnormal_reads2_plus\tnormal_reads2_minus");
		 		if(params.containsKey("output-vcf"))
		 		{
					outValidation.println(vcfHeader);
		 		}
	 		}


	 		BufferedReader normal = new BufferedReader(new FileReader(normalPileupFile));
		    BufferedReader tumor = new BufferedReader(new FileReader(tumorPileupFile));

	    	// If input file not ready, give it a few seconds //
	    	int numNaps = 0;

	    	while(!(normal.ready() && tumor.ready()))
	    	{
	    		try {
			    	Thread.sleep(5000);
			    	numNaps++;

			    	if(numNaps > 100)
			    	{
			    		System.err.println("Input file(s) were not ready for parsing after 100 5-second cycles! Pileup output may be invalid or too slow");
			    		System.exit(10);
			    	}
		    	}
		    	catch(Exception e)
		    	{

		    	}
	    	}


		    if(!(normal.ready() && tumor.ready()))
		    {
		    	System.err.println("ERROR: Input file(s) not ready for parsing! Pileup output may be invalid or too slow.");
		    	System.exit(10);
		    }

		    String lineNormal;
		    String lineTumor;
		    String chromNormal = "";
		    String chromTumor = "";
		    String prevChromNormal = "";
		    String prevChromTumor = "";
		    String refBase = "";
		    int posNormal = 0;
		    int posTumor = 0;


		    DecimalFormat oneDigit = new DecimalFormat("#0.0");
		    DecimalFormat threeDigits = new DecimalFormat("#0.000");


		    // Get first line of Normal //

		    if((lineNormal = normal.readLine()) != null)
		    {
		    	String[] normalContents = lineNormal.split("\t");

		    	if(normalContents.length > 1)
		    	{
			    	chromNormal = normalContents[0];
			    	posNormal = Integer.parseInt(normalContents[1]);
		    	}
		    }

		    // Loop through lines in tumor //

	    	while ((lineTumor = tumor.readLine()) != null)
	    	{
	    		tumorPositions++;
	    		String[] tumorContents = lineTumor.split("\t");

		    	if(tumorContents.length > 1)
		    	{
			    	chromTumor = tumorContents[0];
			    	posTumor = Integer.parseInt(tumorContents[1]);
		    	}

		    	// Parse normal lines until we get the same chromosome //
		    	boolean flagEOF = false;
		    	boolean normalWasReset = false;

		    	//	Advance in normal file if tumor is changed but normal is not, or if tumor is higher //
		    	while(!chromNormal.equals(chromTumor) && !chromTumor.equals(prevChromTumor) && !flagEOF && (chromNormal.equals(prevChromTumor) || inSortOrder(chromNormal, chromTumor)))
		    	{
		    		//System.err.println("Normal (" + chromNormal + ") catching up to " + chromTumor);
		    		// Get next line from normal pileup //
		    		if((lineNormal = normal.readLine()) != null)
		    		{
		    			String[] normalContents = lineNormal.split("\t");

				    	if(normalContents.length > 1)
				    	{
					    	chromNormal = normalContents[0];
					    	posNormal = Integer.parseInt(normalContents[1]);
				    	}
		    		}
		    		else
		    		{
		    			flagEOF = true;
		    		}


		    	}

		    	// If chromosomes match and are non-blank, attempt to get matching positions //
		    	if(chromNormal.equals(chromTumor) && !chromNormal.equals(""))
		    	{
		    		normalWasReset = false;
		    		// Seek to matching Normal Position //

		    		while(chromNormal.equals(chromTumor) && posNormal < posTumor && ((lineNormal = normal.readLine()) != null))
		    		{
		    			String[] normalContents = lineNormal.split("\t");
				    	if(normalContents.length > 1)
				    	{
					    	chromNormal = normalContents[0];
					    	posNormal = Integer.parseInt(normalContents[1]);
				    	}
		    		}

		    		// Seek to matching Tumor Position //

		    		while(chromNormal.equals(chromTumor) && posTumor < posNormal && ((lineTumor = tumor.readLine()) != null))
		    		{
		    			tumorContents = lineTumor.split("\t");
				    	if(tumorContents.length > 1)
				    	{
					    	chromTumor = tumorContents[0];
					    	posTumor = Integer.parseInt(tumorContents[1]);
				    	}
		    		}

		    		// Proceed if normal and tumor positions match //

		    		if(chromNormal.equals(chromTumor) && chromNormal.equals(chromTumor) && posNormal == posTumor)
		    		{
		    			//stats.put("sharedPositions", (stats.get("sharedPositions") + 1));
		    			sharedPositions++;
		    			if(params.containsKey("verbose"))
		    					System.err.println("Comparing calls at " + chromTumor + ":" + posTumor);

		    			refBase = tumorContents[2];
		    			String compareResult = comparePositions(lineNormal, lineTumor, minCoverage, minReads2, minVarFreq, minAvgQual, pValueThreshold, somaticPvalue, minFreqForHom, normalPurity, tumorPurity);

		    			if(compareResult.length() > 0)
		    			{
		    				// Get the alleles to determine type //
		    				String[] compareContents = compareResult.split("\t");
			    			String allele1 = compareContents[0];
			    			String allele2 = compareContents[1];

			    			double strandedness1 = 0.50;
			    			double strandedness2 = 0.50;
			    			double strandednessDiff = 0.00;

			    			if(compareContents.length >= 17)
			    			{
			    				try
			    				{
			    					int tumorReads1plus = Integer.parseInt(compareContents[13]);
			    					int tumorReads1minus = Integer.parseInt(compareContents[14]);
			    					int tumorReads2plus = Integer.parseInt(compareContents[15]);
			    					int tumorReads2minus = Integer.parseInt(compareContents[16]);

			    					if(tumorReads1plus > 0 || tumorReads1minus > 0)
			    					{
			    						strandedness1 = (double) tumorReads1plus / (double) (tumorReads1plus + tumorReads1minus);
			    					}

			    					if(tumorReads2plus > 0 || tumorReads2minus > 0)
			    					{
			    						strandedness2 = (double) tumorReads2plus / (double) (tumorReads2plus + tumorReads2minus);
			    						if(tumorReads1plus > 0 || tumorReads1minus > 0)
			    						{
			    							strandednessDiff = java.lang.Math.abs(strandedness1 - strandedness2);
			    						}
			    					}
			    				}
			    				catch(Exception e)
			    				{
			    					// Exception parsing info from compareResult //
			    				}
			    			}

		    				//stats.put("comparedPositions", (stats.get("comparedPositions") + 1));
			    			comparedPositions++;

		    				if(params.containsKey("verbose") && !compareResult.contains("Reference"))
		    					System.err.println(chromNormal + "\t" + posNormal + "\t" + compareResult);

		    				// If VCF format specified, supply it //

		    				if(params.containsKey("output-vcf"))
		    				{
		    					int normalReads1 = Integer.parseInt(compareContents[2]);
		    					int normalReads2 = Integer.parseInt(compareContents[3]);
		    					String normalFreq = compareContents[4];
		    					String normalCall = compareContents[5];
		    					int tumorReads1 = Integer.parseInt(compareContents[6]);
		    					int tumorReads2 = Integer.parseInt(compareContents[7]);
		    					String tumorFreq = compareContents[8];
		    					String tumorCall = compareContents[9];
		    					String somStatus = compareContents[10];
		    					Double germlineP = Double.parseDouble(compareContents[11]);
		    					Double somaticP = Double.parseDouble(compareContents[12]);

		    					String[] normalContents = lineNormal.split("\t");
		    					//tumorContents = lineTumor.split("\t");
		    					int pileupDepthNormal = Integer.parseInt(normalContents[3]);
		    					int pileupDepthTumor = Integer.parseInt(tumorContents[3]);

		    					int totalDepth = pileupDepthNormal + pileupDepthTumor;

		    					if(allele2.startsWith("+"))
		    					{
		    						// INSERTION //
		    						// Ref = ref base; Var = ref base followed by inserted bases //
		    						String varColumn = allele1 + allele2.replace("+", "");
		    						compareResult = "." + "\t" + allele1 + "\t" + varColumn + "\t" + ".";
		    					}
		    					else if(allele2.startsWith("-"))
		    					{
		    						// DELETION //
		    						// Ref = ref base followed by deleted bases; var = ref base //
		    						String refColumn = allele1 + allele2.replace("-", "");
		    						compareResult = "." + "\t" + refColumn + "\t" + allele1 + "\t" + ".";
		    					}
		    					else
		    					{
			    					compareResult = "." + "\t" + allele1 + "\t" + allele2 + "\t" + ".";
		    					}


		    					// Decide on filter field //
		    					if(params.containsKey("strand-filter") && strandednessDiff > 0.10 && (strandedness2 < 0.10 || strandedness2 > 0.90))
		    					{
		    						compareResult += "\t" + "str10";
		    					}
		    					else if(somStatus.equals("IndelFilter"))
		    					{
		    						compareResult += "\t" + "indelError";
		    					}
		    					else
		    					{
		    						compareResult += "\t" + "PASS";
		    					}

		    					// Determine somatic status id and score //
		    					int ssCode = 0;
		    					double somScore = 0;

		    					if(somStatus.equals("Reference"))
		    					{
		    						// Wildtype //
		    						ssCode = 0;
		    						calledReference++;
		    					}
		    					else if(somStatus.equals("Germline"))
		    					{
		    						// Germline //
		    						ssCode = 1;
		    						calledGermline++;
		    						if(somaticP == 0)
		    						{
		    							somScore = 0;
		    						}
		    						else
		    						{
		    							somScore = 0 - (10 * java.lang.Math.log10(somaticP));
		    						}
		    					}
		    					else if(somStatus.equals("Somatic"))
		    					{
		    						// Somatic //
		    						ssCode = 2;
		    						calledSomatic++;
		    						if(somaticP == 0)
		    						{
		    							somScore = 255;
		    						}
		    						else
		    						{
		    							somScore = 0 - (10 * java.lang.Math.log10(somaticP));
		    						}
		    					}
		    					else if(somStatus.equals("LOH"))
		    					{
		    						// LOH //
		    						ssCode = 3;
		    						calledLOH++;
		    						if(somaticP == 0)
		    						{
		    							somScore = 255;
		    						}
		    						else
		    						{
		    							somScore = 0 - (10 * java.lang.Math.log10(somaticP));
		    						}
		    					}
		    					else
		    					{
		    						// Unknown //
		    						calledUnknown++;
		    						ssCode = 5;
		    					}

		    					// Adjust somatic score //
		    					if(somScore > 255)
		    						somScore = 255;

		    					// Print the info field //

		    					compareResult += "\t" + "DP=" + totalDepth;
		    					if(somStatus.equals("Somatic"))
		    						compareResult += ";SOMATIC";
		    					compareResult += ";" + "SS=" + ssCode;
		    					compareResult += ";" + "SSC=" + (int) somScore;
		    					compareResult += ";" + "GPV=" + pvalueFormat.format(germlineP);
		    					compareResult += ";" + "SPV=" + pvalueFormat.format(somaticP);

		    					// Print the format field //

		    					String tumorDP4 = "";
		    					String normalDP4 = "";

				    			if(compareContents.length >= 17)
				    			{
				    				try
				    				{
				    					tumorDP4 = compareContents[13] + "," + compareContents[14] + "," + compareContents[15] + "," + compareContents[16];
				    					normalDP4 = compareContents[17] + "," + compareContents[18] + "," + compareContents[19] + "," + compareContents[20];
				    				}
				    				catch(Exception e)
				    				{
				    					// Exception parsing info from compareResult //
				    					tumorDP4 = "";
				    					normalDP4 = "";
				    				}
				    			}

				    			if(tumorDP4.length() > 0)
				    				compareResult += "\tGT:GQ:DP:RD:AD:FREQ:DP4";
				    			else
				    				compareResult += "\tGT:GQ:DP:RD:AD:FREQ";

		    					// Determine normal genotype //
		    					String normalGt = ".";
		    					String tumorGt = ".";
		    					if(normalCall.equals(refBase))
		    					{
		    						normalGt = "0/0";
		    					}
		    					else if(VarScan.isHeterozygous(normalCall))
		    					{
		    						normalGt = "0/1";
		    					}
		    					else
		    					{
		    						normalGt = "1/1";
		    					}

		    					if(tumorCall.equals(refBase))
		    					{
		    						tumorGt = "0/0";
		    					}
		    					else if(VarScan.isHeterozygous(tumorCall))
		    					{
		    						tumorGt = "0/1";
		    					}
		    					else
		    					{
		    						tumorGt = "1/1";
		    					}

		    					if(tumorDP4.length() > 0)
		    					{
			    					compareResult += "\t" + normalGt + ":.:" + pileupDepthNormal + ":" + normalReads1 + ":" + normalReads2 + ":" + normalFreq + ":" + normalDP4;
			    					compareResult += "\t" + tumorGt + ":.:" + pileupDepthTumor + ":" + tumorReads1 + ":" + tumorReads2 + ":" + tumorFreq + ":" + tumorDP4;
		    					}
		    					else
		    					{
			    					compareResult += "\t" + normalGt + ":.:" + pileupDepthNormal + ":" + normalReads1 + ":" + normalReads2 + ":" + normalFreq;
			    					compareResult += "\t" + tumorGt + ":.:" + pileupDepthTumor + ":" + tumorReads1 + ":" + tumorReads2 + ":" + tumorFreq;
		    					}
		    				}
		    				// Print to master file for validation //

		    				if(params.containsKey("validation"))
		    				{
		    					outValidation.println(chromNormal + "\t" + posNormal + "\t" + compareResult);
		    				}

		    				if(!params.containsKey("validation") && (compareResult.contains("Reference") || compareResult.contains("SS=0") || compareResult.contains("Filter")))
		    				{
		    					// Don't print reference/indelfilter positions unless doing validation //
		    				}
		    				else if(params.containsKey("strand-filter") && strandednessDiff > 0.10 && (strandedness2 < 0.10 || strandedness2 > 0.90))
		    				{
		    					// If filter is set, ignore variants that are supported largely by one strand //
		    					compareResult = "StrandFilter";
		    				}
		    				else if(allele1.contains("-") || allele1.contains("+") || allele2.contains("-") || allele2.contains("+"))//if(compareResult.contains("INS") || compareResult.contains("DEL"))
		    				{
		    					outIndel.println(chromNormal + "\t" + posNormal + "\t" + compareResult);
		    				}
		    				else
		    				{
			    				outSnp.println(chromNormal + "\t" + posNormal + "\t" + compareResult);
		    				}
		    			}
		    			else
		    			{
//		    				System.err.println("Uncalled" + chromNormal + "\t" + posNormal + "\t" + compareResult);
		    			}

		    			if(compareResult.contains("Reference"))
		    				calledReference++; //stats.put("calledReference", (stats.get("calledReference") + 1));
		    			else if(compareResult.contains("IndelFilter"))
		    				indelFilter++;	//stats.put("indelFilter", (stats.get("indelFilter") + 1));
		    			else if(compareResult.contains("StrandFilter"))
		    				strandFilter++;
		    			else if(compareResult.contains("Germline"))
		    				calledGermline++;	//stats.put("calledGermline", (stats.get("calledGermline") + 1));
		    			else if(compareResult.contains("Somatic"))
		    				calledSomatic++;	//stats.put("calledSomatic", (stats.get("calledSomatic") + 1));
		    			else if(compareResult.contains("LOH"))
		    				calledLOH++;	//stats.put("calledLOH", (stats.get("calledLOH") + 1));
		    			else if(compareResult.contains("Unknown"))
		    				calledUnknown++;	//stats.put("calledUnknown", (stats.get("calledUnknown") + 1));
		    			else if(compareResult.contains("Variant"))
		    				calledVariant++;	//stats.put("calledVariant", (stats.get("calledVariant") + 1));

				    	prevChromNormal = chromNormal;
				    	prevChromTumor = chromTumor;
		    		}
		    		else
		    		{
		    			//System.err.println("Failed to match positions " + chromNormal + " " + posNormal + " to Tumor " + chromTumor + " " + posTumor);
		    		}
		    	}
		    	// If they're in sort order, do nothing so that tumor can catch up //
		    	else if(inSortOrder(chromNormal, chromTumor))
		    	{
		    		System.err.println("Not resetting normal file because " + chromNormal + " < " + chromTumor);
		    	}
		    	// If we reached the end of the normal file but never saw this chromosome, //
		    	// fast-forward until tumor chromosome changes and reset normal file //
		    	else if(flagEOF)
		    	{
		    		flagEOF = false;

		    		while(prevChromTumor.equals(chromTumor) && !flagEOF)
		    		{
		    			if((lineTumor = tumor.readLine()) != null)
		    			{
			    			tumorContents = lineTumor.split("\t");

					    	if(tumorContents.length > 1)
					    	{
						    	chromTumor = tumorContents[0];
						    	posTumor = Integer.parseInt(tumorContents[1]);
					    	}
		    			}
		    			else
		    			{
		    				flagEOF = true;
		    			}
		    		}

		    		// Reset the normal file if we've already passed this chromosome in normal //

		    		if(!flagEOF && !normalWasReset)
		    		{
		    			if(inSortOrder(chromNormal, chromTumor))
		    			{
		    				System.err.println("Not resetting normal file because " + chromNormal + " < " + chromTumor);
		    			}
		    			else
		    			{
		    				System.err.println("Resetting normal file because " + chromNormal + " > " + chromTumor);
				    		normalWasReset = true;
			    			normal.close();
				    		normal = new BufferedReader(new FileReader(normalPileupFile));
		    			}

		    		}
		    	}

	    	}


		    normal.close();
		    tumor.close();



		    outSnp.close();
		    outIndel.close();

		    System.err.println(tumorPositions + " positions in tumor");
		    System.err.println(sharedPositions + " positions shared in normal"); //stats.get("sharedPositions")
		    System.err.println(comparedPositions + " had sufficient coverage for comparison"); //stats.get("comparedPositions")
		    System.err.println(calledReference + " were called Reference"); //stats.get("calledReference")
		    System.err.println(indelFilter + " were mixed SNP-indel calls and filtered");
		    if(params.containsKey("strand-filter"))
		    	System.err.println(strandFilter + " were removed by the strand filter");
		    System.err.println(calledGermline + " were called Germline");
		    System.err.println(calledLOH + " were called LOH");
		    System.err.println(calledSomatic + " were called Somatic");
		    System.err.println(calledUnknown + " were called Unknown");
		    System.err.println(calledVariant + " were called Variant");

		}
		catch (IOException e)
		{
			System.err.println("File Parsing Exception: " + e.getLocalizedMessage());
			e.printStackTrace(System.err);
			System.exit(11);
		}
	}


	/**
	 * Determine if tumor chromosome is before normal chromosome in sort order
	 *
	 * @param	args	Command-line arguments
	 * @return			HashMap of parameter names and their values
	 */
	static boolean inSortOrder(String chrom1, String chrom2)
	{
		String[] testArray = {chrom1, chrom2};
		Arrays.sort(testArray);

		if(testArray[0].equals(chrom1))
			return true;

		return false;
	}

	/**
	 * Parses and verifies any command-line parameters
	 *
	 * @param	args	Command-line arguments
	 * @return			HashMap of parameter names and their values
	 */
	static String comparePositions(String lineNormal, String lineTumor, int minCoverage, int minReads2, double minVarFreq, int minAvgQual, double pValueThreshold, double somaticPvalue, double minFreqForHom, double normalPurity, double tumorPurity)
	{
		try
		{
			 DecimalFormat df = new DecimalFormat("###.##");
			 // Set default parameters //

			 String refBase = "";
			 int normalDepth = 0;
			 String normalBases = "";
			 String normalQualities = "";
			 String normalMapQuals = "";

			 int tumorDepth = 0;
			 String tumorBases = "";
			 String tumorQualities = "";
			 String tumorMapQuals = "";

			 // Parse out normal info //
			 String[] normalContents = lineNormal.split("\t");
			 refBase = normalContents[2].toUpperCase();
			 // Parse out tumor info //
			 String[] tumorContents = lineTumor.split("\t");

			 // Parse out normal from pileup or CNS file //

			 if(normalContents.length >= 6 && normalContents.length <= 7)
			 {
				 normalDepth = Integer.parseInt(normalContents[3]);
				 normalBases = normalContents[4];
				 normalQualities = normalContents[5];
				 normalMapQuals = "";
				 if(normalContents.length == 7)
				 {
					 normalMapQuals = normalContents[6];
				 }

			 }
			 else if(normalContents.length >= 10 && normalContents.length <= 11)
			 {
				 normalDepth = Integer.parseInt(normalContents[7]);
				 normalBases = normalContents[8];
				 normalQualities = normalContents[9];
				 normalMapQuals = "";
				 if(normalContents.length == 11)
				 {
					 normalMapQuals = normalContents[10];
				 }
			 }

			 // Parse out tumor from pileup or CNS file //

			 if(tumorContents.length >= 6 && tumorContents.length <= 7)
			 {
				 tumorDepth = Integer.parseInt(tumorContents[3]);
				 tumorBases = tumorContents[4];
				 tumorQualities = tumorContents[5];
				 tumorMapQuals = "";
				 if(tumorContents.length == 7)
				 {
					 tumorMapQuals = tumorContents[6];
				 }

			 }
			 else if(tumorContents.length >= 10 && tumorContents.length <= 11)
			 {
				 tumorDepth = Integer.parseInt(tumorContents[7]);
				 tumorBases = tumorContents[8];
				 tumorQualities = tumorContents[9];
				 tumorMapQuals = "";
				 if(tumorContents.length == 11)
				 {
					 tumorMapQuals = tumorContents[10];
				 }
			 }


			String somaticStatus = "";
			String allele2 = "";
			double pValue = 1;
			double diffPvalue = 1;

			 if(tumorDepth >= minCoverage && normalDepth >= minCoverage)
			 {
				 // Adjust for tumor purity (i.e., tumor cellularity content of sample) //
				 double tumorMinVarFreq = minVarFreq;

				 // If tumor purity is less than 100%, reduce the minimum variant allele frequency accordingly //
				 if(tumorPurity < 1.00)
				 {
					 tumorMinVarFreq = (minVarFreq * tumorPurity);
				 }

				 HashMap<String, String> readCountsTumor = VarScan.getReadCounts(refBase, tumorBases, tumorQualities, minAvgQual, tumorMapQuals);
	//			 String tumorConsensusLine = callConsensus(refBase, tumorPileup, min_reads2, min_var_freq, min_avg_qual, pValue, purityNormal);
				 String tumorConsensusLine = VarScan.callPosition(refBase, readCountsTumor, "CNS", minReads2, tumorMinVarFreq, minAvgQual, 0.99, minFreqForHom);
				 String[] tumorConsensusContents = tumorConsensusLine.split("\t");
				 String tumorConsensus = tumorConsensusContents[0];

				 if(tumorConsensus.equals("N"))
				 {
					 // No tumor call made, so make no call //
					 return("");
				 }
				 else if(normalDepth >= minCoverage)
				 {
					 // Adjust for normal purity (i.e., tumor contamination of normal in AML) //
					 double normalMinVarFreq = minVarFreq;

					 if(normalPurity < 1.00)
					 {
						 normalMinVarFreq = (minVarFreq / normalPurity);
					 }

					 HashMap<String, String> readCountsNormal = VarScan.getReadCounts(refBase, normalBases, normalQualities, minAvgQual, normalMapQuals);
					 String normalConsensusLine = VarScan.callPosition(refBase, readCountsNormal, "CNS", minReads2, normalMinVarFreq, minAvgQual, 0.99, minFreqForHom); //pValueThreshold, minFreqForHom);

					 String[] normalConsensusContents = normalConsensusLine.split("\t");
					 String normalConsensus = normalConsensusContents[0];

					 if(normalConsensus.equals("N"))
					 {
						 // Make no call at this position //
						 return("");
					 }
					 else
					 {
						 //	Parse out the read counts in tumor //
							int tumorReads1 = Integer.parseInt(tumorConsensusContents[1]);
							int tumorReads2 = Integer.parseInt(tumorConsensusContents[2]);
							int tumorCoverage = tumorReads1 + tumorReads2;
							String tumorAllele2 = VarScan.getVarAllele(refBase, tumorConsensusContents[0]);

							// Parse out strand support in tumor //
							int tumorReads1plus = 0;
							int tumorReads1minus = 0;
							int tumorReads2plus = 0;
							int tumorReads2minus = 0;
							if(tumorConsensusContents.length > 14)
							{
								tumorReads1plus = Integer.parseInt(tumorConsensusContents[11]);
								tumorReads1minus = Integer.parseInt(tumorConsensusContents[12]);
								tumorReads2plus = Integer.parseInt(tumorConsensusContents[13]);
								tumorReads2minus = Integer.parseInt(tumorConsensusContents[14]);
							}

							// Parse out strand support in normal //
							int normalReads1plus = 0;
							int normalReads1minus = 0;
							int normalReads2plus = 0;
							int normalReads2minus = 0;
							if(normalConsensusContents.length > 14)
							{
								normalReads1plus = Integer.parseInt(normalConsensusContents[11]);
								normalReads1minus = Integer.parseInt(normalConsensusContents[12]);
								normalReads2plus = Integer.parseInt(normalConsensusContents[13]);
								normalReads2minus = Integer.parseInt(normalConsensusContents[14]);
							}

							// Parse out the read counts in normal //

							int normalReads1 = Integer.parseInt(normalConsensusContents[1]);
							int normalReads2 = Integer.parseInt(normalConsensusContents[2]);
							int normalCoverage = normalReads1 + normalReads2;
							String normalAllele2 = VarScan.getVarAllele(refBase, normalConsensusContents[0]);


							// Get the Normal Read counts for the tumor variant allele //

							if(!tumorAllele2.equals(refBase)) // normalAllele2.equals(refBase) &&
							{
								allele2 = tumorAllele2;
								if(readCountsNormal.containsKey(tumorAllele2))
								{
									String[] alleleContents = readCountsNormal.get(tumorAllele2).split("\t");
									normalReads2 = Integer.parseInt(alleleContents[0]);
									normalCoverage = normalReads1 + normalReads2;
								}
							}
							else if(!normalAllele2.equals(refBase))
							{
								allele2 = normalAllele2;
							}
							else
							{
								// Neither consensus contained a variant allele, so get most-observed tumor variant //
								if(tumorConsensusContents.length > 15)
								{
									allele2 = tumorConsensusContents[15];
								}
								else if(tumorConsensusContents.length == 10)
								{
									allele2 = tumorConsensusContents[9];
								}
								else if(normalConsensusContents.length > 15)
								{
									allele2 = normalConsensusContents[15];
								}
								else if(normalConsensusContents.length == 10)
								{
									allele2 = normalConsensusContents[9];
								}
							}


							double normalFreq = (double) normalReads2 / (double) normalCoverage;
							double tumorFreq = (double) tumorReads2 / (double) tumorCoverage;

							// Calculate the frequency difference //
							double freqDiff = tumorFreq - normalFreq;

							 // P-value of significant difference //
							diffPvalue = VarScan.getSignificance(normalReads1, normalReads2, tumorReads1, tumorReads2);

							// Format allele frequencies for printing //
							String normalFreqPrint = df.format(normalFreq * 100) + "%";
							String tumorFreqPrint = df.format(tumorFreq * 100) + "%";

							 // If Normal matches Tumor it's either reference or Germline //

							 if(normalConsensus.equals(tumorConsensus) && (normalConsensus.equals(refBase) || diffPvalue > somaticPvalue))
							 {
								 // CASE 0: Normal and Tumor Match //

								 if(normalConsensus.equals(refBase))
								 {
									 somaticStatus = "Reference";
								 }
								 else
								 {
									 // Recalculate p-value //
									 int totalReads1 = normalReads1 + tumorReads1;
									 int totalReads2 = normalReads2 + tumorReads2;
									 int totalCoverage = totalReads1 + totalReads2;
									 // P-value of Germline variant //
									 pValue = VarScan.getSignificance(totalCoverage, 0, totalReads1, totalReads2);
									 if(pValue <= somaticPvalue) //Changed from if(pValue <= pValueThreshold) 11-jun-2012
									 {
										 somaticStatus = "Germline";
										 allele2 = tumorAllele2;
									 }
									 else
									 {
										 somaticStatus = "Reference";
										 allele2 = refBase;
									 }
								 }
							 }

							 // If Normal does NOT match Tumor it could be Somatic, LOH, or Unknown //

							 else
							 {
								 if(normalConsensus.equals(tumorConsensus))
								 {
									 // Genotype calls match, but the difference must have been significant. //
									 // Let's try harder to call a variant for tumor here //
				//					 tumorConsensusLine = VarScan.callPosition(refBase, readCountsTumor, "CNS", 1, 0.00, 0, 1.00, minFreqForHom);
				//					 tumorConsensusContents = tumorConsensusLine.split("\t");
				//					 tumorConsensus = tumorConsensusContents[0];
				//					 tumorAllele2 = VarScan.getVarAllele(refBase, tumorConsensusContents[0]);
				//					 System.err.println("Got a new consensus: " + tumorConsensus + " from " + tumorAllele2);
								 }
								 // CASE 1: Indel-associated SNP Filter //

								 if(tumorConsensus.contains("/") && !normalConsensus.contains("/") && !normalConsensus.equals(refBase))
								 	somaticStatus = "IndelFilter";
								 else if(normalConsensus.contains("/") && !tumorConsensus.contains("/") && !tumorConsensus.equals(refBase))
								 	somaticStatus = "IndelFilter";
								 else
								 {
				//					 CASE 2: Somatic indel or SNP events, where difference in read counts is significant or else coverage is low //

										if(diffPvalue <= somaticPvalue || normalFreq == 0.00) // || tumorCoverage < 30 || normalCoverage < 15)
										{
											// CASE 2A: Perfect Somatic Het //
											if(normalConsensus.equals(refBase) && VarScan.isHeterozygous(tumorConsensus) && tumorFreq > normalFreq)
											{
												somaticStatus = "Somatic";
												allele2 = tumorAllele2;
											}
											// CASE 2B: Somatic Homozygous //
											else if(normalConsensus.equals(refBase) && VarScan.isHomozygous(tumorConsensus))
											{
												somaticStatus = "Somatic";
												allele2 = tumorAllele2;
											}
											// CASE 2C: LOH of variant allele //
											else if(tumorConsensus.equals(refBase) && VarScan.isHeterozygous(normalConsensus))
											{
												somaticStatus = "LOH";
												allele2 = normalAllele2;
											}
											// CASE 2D: LOH of reference allele //
											else if(VarScan.isHeterozygous(normalConsensus) && VarScan.isHomozygous(tumorConsensus))
											{
												somaticStatus = "LOH";
												allele2 = tumorAllele2;
											}
											// CASE 2E: Variant alleles match but difference significant //
											else if(tumorAllele2.equals(normalAllele2))
											{
												if(normalFreq > minVarFreq)
												{
													somaticStatus = "Germline";
												}
												else if(freqDiff >= 0.30 && tumorFreq > normalFreq)
												{
													somaticStatus = "Somatic";
												}
												else if(freqDiff <= -0.30 && tumorFreq < normalFreq)
												{
													somaticStatus = "LOH";
												}
												else// if(freqDiff < 0.50)
												{
													somaticStatus = "Germline"; // Should this be GOH? //
													// Recalculate p-value //
													int totalReads1 = normalReads1 + tumorReads1;
													int totalReads2 = normalReads2 + tumorReads2;
													int totalCoverage = totalReads1 + totalReads2;
													pValue = VarScan.getSignificance(totalCoverage, 0, totalReads1, totalReads2);
												}
				//								else
				//								{
				//									somaticStatus = "LOH";
				//								}
												allele2 = tumorAllele2;
											}
											// CASE 2F: Variant alleles don't match but tumor het and higher freq = normal FalsePos //
											else if(tumorFreq > normalFreq && VarScan.isHeterozygous(normalConsensus) && VarScan.isHeterozygous(tumorConsensus))
											{
												normalConsensus = refBase;
												somaticStatus = "Somatic";
												allele2 = tumorAllele2;
											}
											// CASE 2G: Unknown Somatic Change, e.g. reverse-LOH (GOH) //
											else
											{
												somaticStatus = "Unknown";
												if(tumorAllele2.equals(refBase))
													allele2 = normalAllele2;
												else
													allele2 = tumorAllele2;
											}

										}
										else
										{
											// CASE 3: Difference not significant //

											// CASE 3A: One sample het, one sample hom = Germline //
											if(tumorAllele2.equals(normalAllele2))
											{
												// Recalculate p-value //
												int totalReads1 = normalReads1 + tumorReads1;
												int totalReads2 = normalReads2 + tumorReads2;
												int totalCoverage = totalReads1 + totalReads2;
												pValue = VarScan.getSignificance(totalCoverage, 0, totalReads1, totalReads2);
												if(pValue <= pValueThreshold)
												{
													somaticStatus = "Germline";
													allele2 = tumorAllele2;
												}
												else
												{
													somaticStatus = "Reference";
													allele2 = refBase;
												}
											}

											// CASE 3B: Probable false positive in tumor//
											else if(normalConsensus.equals(refBase))
											{
												somaticStatus = "Reference";
												allele2 = tumorAllele2;
											}

											// CASE 3C: Probable false positive in numor//
											else if(tumorConsensus.equals(refBase))
											{
												somaticStatus = "Reference";
												allele2 = normalAllele2;
											}
											else
											{
												somaticStatus = "Unknown";
												allele2 = normalAllele2 + "/" + tumorAllele2;
											}

										}

									 }

								 }


							 // Compile the report //

							 String resultLine = refBase + "\t" + allele2 + "\t";
							 resultLine += normalReads1 + "\t" + normalReads2 + "\t" + normalFreqPrint + "\t" + normalConsensus + "\t";
							 resultLine += tumorReads1 + "\t" + tumorReads2 + "\t" + tumorFreqPrint + "\t" + tumorConsensus + "\t";
							 resultLine += somaticStatus + "\t" + pValue + "\t" + diffPvalue + "\t";
							 resultLine += tumorReads1plus + "\t" + tumorReads1minus + "\t";
							 resultLine += tumorReads2plus + "\t" + tumorReads2minus + "\t";
							 resultLine += normalReads1plus + "\t" + normalReads1minus + "\t";
							 resultLine += normalReads2plus + "\t" + normalReads2minus;
							 return(resultLine);

					 }


				 }
				 else
				 {
					 return(""); // Normal did not meet coverage
				 }

			}
			else
			{
				 return(""); // Tumor did not meet coverage
			}
		}
		catch(Exception e)
		{
			System.err.println("Warning:");
			e.printStackTrace(System.err);
		}
		 return("");	// No call
	}



	/**
	 * Determines the sort order for chromosomes
	 *
	 * @param	args	Command-line arguments
	 * @return			HashMap of parameter names and their values
	 */
	static Boolean chromSorted(String chrom1, String chrom2)
	{
		Boolean answer = false;

		chrom1.replace("X", "23");
		chrom1.replace("Y", "24");
		chrom1.replace("M", "25");

		chrom2.replace("X", "23");
		chrom2.replace("Y", "24");
		chrom2.replace("M", "25");

		String[] unsorted = {chrom1, chrom2};
		String[] sorted = {chrom1, chrom2};
		Arrays.sort(sorted);
		System.err.println("Sorted order is " + sorted[0] + " " + sorted[1]);
		try{
			if(sorted[0].equals(unsorted[0]))
			{
				answer = true;
			}
		}
		catch(Exception e)
		{

		}

		return(answer);
	}

}
