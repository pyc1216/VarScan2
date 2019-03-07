/**
 * @(#)CallPileup.java
 *
 * Copyright (c) 2009-2010 Daniel C. Koboldt and Washington University in St. Louis
 *
 * COPYRIGHT
 */

package net.sf.varscan;

import java.io.BufferedReader;
import java.util.HashMap;

/**
 * A class for calling variants or consensus bases from a pileup file
 *
 * @version	2.3
 *
 * @author Daniel C. Koboldt <dkoboldt@genome.wustl.edu>
 *
 */
public class CallPileup {

	public CallPileup(String[] args, String callType)
	{
		// Define the usage message //
		String usage = "USAGE: java -jar VarScan.jar pileup2cns [pileup file] OPTIONS\n" +
		"\tpileup file - The SAMtools pileup file\n" +
		"\n" +
		"\tOPTIONS:\n" +
		"\t--min-coverage\tMinimum read depth at a position to make a call [8]\n" +
		"\t--min-reads2\tMinimum supporting reads at a position to call variants [2]\n" +
		"\t--min-avg-qual\tMinimum base quality at a position to count a read [15]\n" +
		"\t--min-var-freq\tMinimum variant allele frequency threshold [0.01]\n" +
		"\t--min-freq-for-hom\tMinimum frequency to call homozygote [0.75]\n" +
		"\t--p-value\tDefault p-value threshold for calling variants [99e-02]\n" +
		"\t--variants\tReport only variant (SNP/indel) positions [0]";

		// Set parameter defaults //

		HashMap<String, String> params = VarScan.getParams(args);

		int minCoverage = 8;
		int minReads2 = 2;
		int minAvgQual = 15;
		double minVarFreq = 0.01;
		double minFreqForHom = 0.75;
		double pValueThreshold = 0.99;

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

//	    HashMap<String, Integer> stats = new HashMap<String, Integer>();
//	    stats.put("numBases", 0);
//	    stats.put("numCovered", 0);
//	    stats.put("numCalled", 0);
//	    stats.put("calledRef", 0);
//	    stats.put("calledIndel", 0);
//	    stats.put("calledSNP", 0);
//	    stats.put("numParsingExceptions", 0);

		long numBases = 0;
		long numCovered = 0;
		long numCalled = 0;
		long calledRef = 0;
		long calledIndel = 0;
		long calledSNP = 0;
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

	    	if(in != null && in.ready())
	    	{
	    		// Print a file header //
	    		if(!params.containsKey("no-headers"))
	    			System.out.println("Chrom\tPosition\tRef\tCons\tReads1\tReads2\tVarFreq\tStrands1\tStrands2\tQual1\tQual2\tPvalue\tMapQual1\tMapQual2\tReads1Plus\tReads1Minus\tReads2Plus\tReads2Minus\tVarAllele");

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

	    				// Verify expected pileup format //

	    				if(lineContents.length > 5 && lineContents[0].length() > 0 && lineContents[1].length() > 0 && lineContents[2].length() > 0 && lineContents[3].length() > 0)
	    				{
	    					String refName = "";
	    					String position = "";
	    					String refBase = "";
	    					int readDepth = 0;
	    					String readBases = "";
	    					String readQualities = "";
	    					String mapQualities = "";

	    					// Pileup Files have 6-7 columns //
	    					if(lineContents.length <= 7)
	    					{
		    					refName = lineContents[0];
			    	        	position = lineContents[1];
			    	        	refBase = lineContents[2].toUpperCase();
			    	        	readDepth = Integer.parseInt(lineContents[3]);
			    	        	readBases = lineContents[4];
			    	        	readQualities = lineContents[5];
			    	        	mapQualities = "";
			    	        	if(lineContents.length > 6)			// Get Map Qualities if available //
			    	        		mapQualities = lineContents[6];
	    					}
	    					// Pileup lines in CNS files have 10-11 columns
	    					else if (lineContents.length >= 10 && lineContents.length <= 11)
	    					{
		    					refName = lineContents[0];
			    	        	position = lineContents[1];
			    	        	refBase = lineContents[2].toUpperCase();
			    	        	readDepth = Integer.parseInt(lineContents[7]);
			    	        	readBases = lineContents[8];
			    	        	readQualities = lineContents[9];
			    	        	mapQualities = "";
			    	        	if(lineContents.length > 10)			// Get Map Qualities if available //
			    	        		mapQualities = lineContents[10];
	    					}
	    					// Indel calls in CNS files have 15-16 columns
	    					else if(lineContents.length >= 15 && lineContents.length <= 16)
	    					{
	    						// Ignore these //
	    					}


		    	        	if(readDepth >= minCoverage && VarScan.qualityDepth(readQualities, minAvgQual) >= minCoverage)
		    	        	{
		    	        		numCovered++;

		    	        		// Build brief pileup string //

		    	        		HashMap<String, String> readCounts = VarScan.getReadCounts(refBase, readBases, readQualities, minAvgQual, mapQualities);

		    	        		String positionCall = VarScan.callPosition(refBase, readCounts, callType, minReads2, minVarFreq, minAvgQual, pValueThreshold, minFreqForHom);

		    	        		if(positionCall.length() > 0)
		    	        		{
	    	        				numCalled++;
		    	        			String[] callLines = positionCall.split("\n");

		    	        			// Go thru each line in resulting call list //
		    	        			for(int lineCounter = 0; lineCounter < callLines.length; lineCounter++)
		    	        			{
		    	        				// Determine type of call that was made //
		    	        				String[] callContents = callLines[lineCounter].split("\t");
		    	        				String consBase = callContents[0];

		    	        				if(consBase.equals(refBase))
		    	        				{
		    	        					calledRef++;
		    	        				}
		    	        				else if(consBase.length() > 1)
		    	        				{
		    	        					calledIndel++;
		    	        				}
		    	        				else
		    	        				{
		    	        					calledSNP++;
		    	        				}

				    	        		// Print some results //
		    	        				if(params.containsKey("variants") && (consBase.equals(refBase) || consBase.equals("N")))
		    	        				{
		    	        					// Don't print ref base if only printing variants //
		    	        				}
		    	        				else
		    	        				{
		    	        					System.out.println(refName + "\t" + position + "\t" + refBase + "\t" + callLines[lineCounter]);
		    	        				}
		    	        			}
		    	        		}


		    	        	}
		    	        	else
		    	        	{
		    	        		// Either raw depth or quality depth did not meet minimum //
		    	        		if(readDepth >= minCoverage)
		    	        		{
		    	        			// Raw depth was enough, but quality depth was not //
		    	        			 //System.err.println("Raw depth = " + readDepth + " Qual depth = " + VarScan.qualityDepth(readQualities, minAvgQual) + " < " + minCoverage);
		    	        		}
		    	        	}
	    				}
	    				else
	    				{
	    					System.err.println("Error: Invalid format for pileup at line " + numBases + "\n" + line + "\n");
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

				in.close();

				System.err.println(numBases + " bases in pileup file");
				System.err.println(numCovered + " met minimum coverage of " + minCoverage + "x");

				if(callType.equals("SNP"))
					System.err.println(calledSNP + " SNPs predicted");
				else if(callType.equals("INDEL"))
					System.err.println(calledIndel + " indels predicted");
				else //CNS //
				{
					System.err.println(numCalled + " positions were called");
					System.err.println(calledRef + " called Reference");
					System.err.println(calledSNP + " called SNP");
					System.err.println(calledIndel + " called indel");
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
}
