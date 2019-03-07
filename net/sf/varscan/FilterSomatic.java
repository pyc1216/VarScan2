/**
 * @(#)FilterSomatic.java
 *
 * Copyright (c) 2009-2010 Daniel C. Koboldt and Washington University in St. Louis
 *
 * COPYRIGHT
 */

package net.sf.varscan;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.HashMap;

/**
 * A class for filtering VarScan variant predictions
 *
 * @version	2.3
 *
 * @author Daniel C. Koboldt <dkoboldt@genome.wustl.edu>
 *
 */
public class FilterSomatic {

	public FilterSomatic(String[] args)
	{
		//		 Define the usage message //
		String usage = "USAGE: java -jar VarScan.jar filter [variant file] OPTIONS\n" +
		"\tvariant file - A file of SNPs or indels\n" +
		"\n" +
		"\tOPTIONS:\n" +
		"\t--min-coverage\tMinimum read depth at a position to make a call [10]\n" +
		"\t--min-reads2\tMinimum supporting reads at a position to call variants [4]\n" +
		"\t--min-strands2\tMinimum # of strands on which variant observed (1 or 2) [1]\n" +
		"\t--min-var-freq\tMinimum variant allele frequency threshold [0.20]\n" +
		"\t--p-value\tDefault p-value threshold for calling variants [5e-02]\n" +
		"\t--indel-file\tFile of indels for filtering nearby SNPs\n" +
		"\t--output-file\tOptional output file for filtered variants";

		// Set parameter defaults //
//		Set parameter defaults //

		int minCoverage = 	8;
		int minReads2 = 	4;
		int minStrands2 = 	1;
		int minAvgQual = 	20;
		double minVarFreq = 0.20;
		double pValueThreshold = 0.05; //1.0e-04;

		int windowSize = 10;
		int windowSNPs = 3;
		int indelMargin = 3;
		String outFileName = "";

		HashMap<String, Boolean> indelPositions = new HashMap<String, Boolean>();

		// Adjust parameters based on user input //

		HashMap<String, String> params = VarScan.getParams(args);

		try
		{
			if(params.containsKey("output-file"))
				outFileName = params.get("output-file");

			if(params.containsKey("window-size"))
				 windowSize = Integer.parseInt(params.get("window-size"));

			if(params.containsKey("window-snps"))
				 windowSNPs = Integer.parseInt(params.get("window-snps"));

			if(params.containsKey("indel-margin"))
				 indelMargin = Integer.parseInt(params.get("indel-margin"));

			if(params.containsKey("indel-file"))
			{
				indelPositions = loadIndels(params.get("indel-file"));
			}

			if(params.containsKey("min-coverage"))
				 minCoverage = Integer.parseInt(params.get("min-coverage"));

			if(params.containsKey("min-reads2"))
				 minReads2 = Integer.parseInt(params.get("min-reads2"));

			if(params.containsKey("min-strands2"))
				 minStrands2 = Integer.parseInt(params.get("min-strands2"));

			if(params.containsKey("min-var-freq"))
				 minVarFreq = Double.parseDouble(params.get("min-var-freq"));

			if(params.containsKey("min-avg-qual"))
				 minAvgQual = Integer.parseInt(params.get("min-avg-qual"));

			if(params.containsKey("p-value"))
				 pValueThreshold = Double.parseDouble(params.get("p-value"));


			System.err.println("Window size:\t" + windowSize);
			System.err.println("Window SNPs:\t" + windowSNPs);
			System.err.println("Indel margin:\t" + indelMargin);

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

		// Identify SNP clusters //

		HashMap<String, Boolean> clusterSNPs = findSNPclusters(windowSize, windowSNPs, args);


	    // Define two-decimal-place format and statistics hash //

	    HashMap<String, Integer> stats = new HashMap<String, Integer>();
	    stats.put("numVariants", 0);
	    stats.put("numFailCoverage", 0);
	    stats.put("numFailVarFreq", 0);
	    stats.put("numFailPvalue", 0);
	    stats.put("numFailReads2", 0);
	    stats.put("numNearIndel", 0);
	    stats.put("numSNPcluster", 0);
	    stats.put("numPassFilter", 0);
	    stats.put("numParsingExceptions", 0);

		// Parse piped input or user-provided pileup file //

	    try
	    {
			// Declare output file //
			PrintStream outFile = null;
			if(params.containsKey("output-file"))
				outFile = new PrintStream( new FileOutputStream(outFileName) );

	    	// Declare file-parsing variables //

	    	BufferedReader in = VarScan.getInfile(args);

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

	    	String line;
	    	int lineCounter = 0;
	    	boolean isVCF = false;

	    	// Proceed if input stream is ready //

	    	if(in != null && in.ready())
	    	{
	    		// Parse the infile line by line //

	    		while ((line = in.readLine()) != null)
	    		{
	    			lineCounter++;

	    			try
	    			{
		    			String[] lineContents = line.split("\t");
		    			String chrom = lineContents[0];

	    				if(line.startsWith("#"))
	    					isVCF = true;

		    			if(chrom.equals("Chrom") || chrom.equals("chrom") || line.startsWith("#"))
		    			{

		    				// Print header //
	    					if(params.containsKey("output-file"))
	    						outFile.println(line);
		    			}
		    			else
		    			{
		    				stats.put("numVariants", (stats.get("numVariants") + 1));
		    				int position = Integer.parseInt(lineContents[1]);
			    			String thisKey = chrom + "\t" + position;

		    				int normalReads1 = 0;
		    				int normalReads2 = 0;
		    				int tumorReads1 = 0;
		    				int tumorReads2 = 0;
		    				String somaticStatus = "";
		    				double somaticPvalue = 0.00;

	    					if(isVCF)
	    					{
	    						String info = lineContents[7];
	    						String normal = lineContents[9];
	    						String tumor = lineContents[10];

	    						String[] infoContents = info.split(";");
	    						for(int colCounter = 0; colCounter < infoContents.length; colCounter++)
	    						{
	    							String element = infoContents[colCounter];
	    							String[] elementContents = element.split("=");
	    							if(elementContents[0].equals("SS"))
	    								somaticStatus = elementContents[1];
	    							else if(elementContents[0].equals("GPV") && somaticStatus.equals("1"))
	    								somaticPvalue = Double.parseDouble(elementContents[1]);
	    							else if(elementContents[0].equals("SPV") && !somaticStatus.equals("1"))
	    								somaticPvalue = Double.parseDouble(elementContents[1]);
	    						}

	    						String[] normalContents = normal.split(":");
	    						normalReads1 = Integer.parseInt(normalContents[3]);
	    						normalReads2 = Integer.parseInt(normalContents[4]);

	    						String[] tumorContents = tumor.split(":");
	    						tumorReads1 = Integer.parseInt(tumorContents[3]);
	    						tumorReads2 = Integer.parseInt(tumorContents[4]);
	    					}
	    					else
	    					{
			    				normalReads1 = Integer.parseInt(lineContents[4]);
			    				normalReads2 = Integer.parseInt(lineContents[5]);
			    				tumorReads1 = Integer.parseInt(lineContents[8]);
			    				tumorReads2 = Integer.parseInt(lineContents[9]);

			    				somaticStatus = lineContents[12];
			    				somaticPvalue = Double.parseDouble(lineContents[14]);
	    					}

	    					// Proceed //

		    				double normalFreq = 0;
		    				double tumorFreq = 0;

		    				int normalCoverage = normalReads1 + normalReads2;
		    				int tumorCoverage = tumorReads1 + tumorReads2;

		    				if(normalReads1 > 0 || normalReads2 > 0)
		    					normalFreq = (double) normalReads2 / (double) (normalReads1 + normalReads2);

		    				if(tumorReads1 > 0 || tumorReads2 > 0)
		    					tumorFreq = (double) tumorReads2 / (double) (tumorReads1 + tumorReads2);


		    				if(normalReads1 > 0 || normalReads2 > 0)
		    					normalFreq = (double) normalReads2 / (double) (normalReads1 + normalReads2);

		    				if(tumorReads1 > 0 || tumorReads2 > 0)
		    					tumorFreq = (double) tumorReads2 / (double) (tumorReads1 + tumorReads2);

		    				boolean filterFlag = false;
		    				boolean indelFlag = false;
		    				boolean clusterFlag = false;

		    				if(normalCoverage < minCoverage || tumorCoverage < minCoverage)
		    				{
		    					// Fail due to coverage //
		    					filterFlag = true;
		    					stats.put("numFailCoverage", (stats.get("numFailCoverage") + 1));
		    				}
		    				else if(tumorReads2 < minReads2)
		    				{
		    					filterFlag = true;
		    					stats.put("numFailReads2", (stats.get("numFailReads2") + 1));
		    				}
		    				else if(tumorFreq < minVarFreq)
		    				{
		    					filterFlag = true;
		    					stats.put("numFailVarFreq", (stats.get("numFailVarFreq") + 1));
		    				}
		    				else if(somaticPvalue > pValueThreshold)
		    				{
		    					filterFlag = true;
		    					stats.put("numFailPvalue", (stats.get("numFailPvalue") + 1));
		    				}

		    				if(clusterSNPs.containsKey(thisKey))
		    				{
		    					clusterFlag = true;
		    				}
		    				// If indel file was provided, check position-1 to position+1 //
		    				else if(params.containsKey("indel-file"))
		    				{
		    					for(int thisPosition = position - indelMargin; thisPosition <= position + indelMargin; thisPosition++)
		    					{
		    						String key = chrom + "\t" + thisPosition;
		    						if(indelPositions.containsKey(key))
		    						{
		    							indelFlag = true;
		    						}
		    					}
		    				}

		    				if(filterFlag)
		    				{
		    					// Already counted //
		    				}
		    				else if(indelFlag)
		    				{
		    					stats.put("numNearIndel", (stats.get("numNearIndel") + 1));
		    				}
		    				else if(clusterFlag)
		    				{
		    					// Remove cluster SNPs //
		    					stats.put("numSNPcluster", (stats.get("numSNPcluster") + 1));
		    				}
		    				else
		    				{
		    					if(params.containsKey("output-file"))
		    						outFile.println(line);

		    					stats.put("numPassFilter", (stats.get("numPassFilter") + 1));
		    				}


		    			}
	    			}
	    			catch(Exception e)
	    			{
	    				if(lineCounter == 1)
	    				{
	    					if(params.containsKey("output-file"))
	    						outFile.println(line);
	    				}
	    				else
	    				{
		    				System.err.println("Parsing Exception on line:\n" + line + "\n" + e.getLocalizedMessage());
		    				stats.put("numParsingExceptions", (stats.get("numParsingExceptions") + 1));
		    				if(stats.get("numParsingExceptions") >= 5)
		    				{
		    					System.err.println("Too many parsing exceptions encountered; exiting");
		    					return;
		    				}
	    				}

	    			}

	    		}

	    		// Report summary of results //

	    		System.err.println(stats.get("numVariants") + " variants in input stream");
	    		System.err.println(stats.get("numFailCoverage") + " failed to meet coverage requirement");
	    		System.err.println(stats.get("numFailReads2") + " failed to meet reads2 requirement");
	    		System.err.println(stats.get("numFailVarFreq") + " failed to meet varfreq requirement");
	    		System.err.println(stats.get("numFailPvalue") + " failed to meet p-value requirement");
	    		System.err.println(stats.get("numSNPcluster") + " in SNP clusters were removed");
	    		System.err.println(stats.get("numNearIndel") + " were removed near indels");
	    		System.err.println(stats.get("numPassFilter") + " passed filters");

	    		in.close();
	    	}
	    	else
	    	{
	    		System.err.println("Input file not found!");
	    		System.err.println(usage);
	    	}
	    }
	    catch(Exception e)
	    {
	    	System.err.println("Error Parsing Input File: " + e.getLocalizedMessage());
	    	e.printStackTrace(System.err);
	    	return;
	    }
	}


	/**
	 * Loads indels to be used for filtering
	 *
	 * @param	filename	Path to file of indels
	 * @return	indels		HashMap of indel positions (chrom\tposition)
	 */
	static HashMap<String, Boolean> loadIndels(String filename)
	{
    	HashMap<String, Boolean> indels = new HashMap<String, Boolean>();

	    try
	    {
	    	// Declare file-parsing variables //

	    	String line;

    		File infile = new File(filename);
    		if(infile.exists())
    		{
    			BufferedReader in = new BufferedReader(new FileReader(infile));

    			if(in.ready())
    			{
    				while ((line = in.readLine()) != null)
    				{
    					String[] lineContents = line.split("\t");
    					String chrom = lineContents[0];
    					if(chrom.equals("Chrom") || line.startsWith("#"))
    					{
    						// Ignore headers //
    					}
    					else
    					{
	    					String position = lineContents[1];
	    					String indelKey = chrom + "\t" + position;
	    					indels.put(indelKey, Boolean.TRUE);
    					}

    				}
    			}
    			else
    			{
    				System.err.println("Unable to open indels file for reading");
    			}

    			in.close();
    		}
	    }
	    catch(Exception e)
	    {
	    	System.err.println("Error Parsing Indel File: " + e.getLocalizedMessage());
	    }

	    return(indels);
	}


	/**
	 * Loads indels to be used for filtering
	 *
	 * @param	filename	Path to file of indels
	 * @return	indels		HashMap of indel positions (chrom\tposition)
	 */
	static HashMap<String, Boolean> findSNPclusters(int windowSize, int windowSNPs, String[] args)
	{
    	HashMap<String, Boolean> snps = new HashMap<String, Boolean>();
    	HashMap<String, Boolean> clusterSNPs = new HashMap<String, Boolean>();

    	// Declare file-parsing variables //

    	BufferedReader in = VarScan.getInfile(args);
    	String line;
    	int lineCounter = 0;

    	// Proceed if input stream is ready //
    	try
    	{
	    	if(in != null && in.ready())
	    	{
	    		// Parse the infile line by line //

	    		while ((line = in.readLine()) != null)
	    		{
	    			lineCounter++;

	    			try
	    			{
		    			String[] lineContents = line.split("\t");
		    			String chrom = lineContents[0];


		    			if(chrom.equals("Chrom") || line.startsWith("#"))
		    			{
		    				// Ignore headers //
		    			}
		    			else
		    			{
			    			int position = Integer.parseInt(lineContents[1]);
	    					String snpKey = chrom + "\t" + position;
	    					snps.put(snpKey, Boolean.TRUE);
		    			}
	    			}
	    			catch(Exception e)
	    			{
	    				// Ignore parsing exceptions for now
	    			}
	    		}

	    		in.close();
	    	}
    	}
    	catch(Exception e)
    	{
    		// File parsing exception //
    		System.err.println("Error loading SNP positions:" + e.getMessage());
    	}

    	int numClusterSNPs = 0;
    	// Go through each position in SNP keys //

		String[] snpKeys = (String[]) snps.keySet().toArray(new String[0]);
		Arrays.sort(snpKeys);

		for(String snpPosition : snpKeys)
		{
			try
			{
				String[] snpContents = snpPosition.split("\t");
				String chrom = snpContents[0];
    			int position = Integer.parseInt(snpContents[1]);

    			// Search window size in each direction //
    			int numSNPsInWindow = 1;

    			// Check downstream window //
    			for(int thisPosition = position + 1; thisPosition <= position + windowSize; thisPosition++)
    			{
    				String thisKey = chrom + "\t" + thisPosition;
    				if(snps.containsKey(thisKey))
    				{
    					numSNPsInWindow++;
    				}
    			}

    			// If we have a cluster, mark this position as well as any downstream //

    			if(numSNPsInWindow >= windowSNPs)
    			{
    				clusterSNPs.put(snpPosition, true);

        			for(int thisPosition = position + 1; thisPosition <= position + windowSize; thisPosition++)
        			{
        				String thisKey = chrom + "\t" + thisPosition;
        				if(snps.containsKey(thisKey))
        				{
        					clusterSNPs.put(thisKey, true);
        					numClusterSNPs++;
        				}
        			}
    			}



			}
			catch(Exception e)
			{

			}
		}


		System.err.println(numClusterSNPs + " cluster SNPs identified");
	    return(clusterSNPs);
	}
}