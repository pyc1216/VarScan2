/**
 * @(#)FilterVariants.java
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
import java.util.HashMap;

/**
 * A class for filtering VarScan variant predictions
 *
 * @version	2.3
 *
 * @author Daniel C. Koboldt <dkoboldt@genome.wustl.edu>
 *
 */
public class FilterVariants {

	public FilterVariants(String[] args)
	{
		//		 Define the usage message //
		String usage = "USAGE: java -jar VarScan.jar filter [variant file] OPTIONS\n" +
		"\tvariant file - A file of SNPs or indels\n" +
		"\n" +
		"\tOPTIONS:\n" +
		"\t--min-coverage\tMinimum read depth at a position to make a call [10]\n" +
		"\t--min-reads2\tMinimum supporting reads at a position to call variants [2]\n" +
		"\t--min-strands2\tMinimum # of strands on which variant observed (1 or 2) [1]\n" +
		"\t--min-avg-qual\tMinimum average base quality for variant-supporting reads [20]\n" +
		"\t--min-var-freq\tMinimum variant allele frequency threshold [0.20]\n" +
		"\t--p-value\tDefault p-value threshold for calling variants [1e-01]\n" +
		"\t--indel-file\tFile of indels for filtering nearby SNPs\n" +
		"\t--output-file\tFile to contain variants passing filters\n";

		// Set parameter defaults //

		int minCoverage = 	10;
		int minReads2 = 	2;
		int minStrands2 = 	1;
		int minAvgQual = 	15;
		double minVarFreq = 0.20;
		double pValueThreshold = 0.10;
		HashMap<String, Boolean> indelPositions = new HashMap<String, Boolean>();
		String outFileName = "";
		String notFileName = "";

		// Adjust parameters based on user input //

		HashMap<String, String> params = VarScan.getParams(args);

		try
		{
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

			if(params.containsKey("indel-file"))
			{
				indelPositions = loadIndels(params.get("indel-file"));
			}

			if(params.containsKey("output-file"))
				outFileName = params.get("output-file");

			if(params.containsKey("not-file"))
				notFileName = params.get("not-file");

			System.err.println("Min coverage:\t" + minCoverage);
			System.err.println("Min reads2:\t" + minReads2);
			System.err.println("Min strands2:\t" + minStrands2);
			System.err.println("Min var freq:\t" + minVarFreq);
			System.err.println("Min avg qual:\t" + minAvgQual);
			System.err.println("P-value thresh:\t" + pValueThreshold);

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

		// Get the input file and parse it //

	    // Define two-decimal-place format and statistics hash //

	    HashMap<String, Integer> stats = new HashMap<String, Integer>();
	    stats.put("numVariants", 0);
	    stats.put("numNearIndel", 0);
	    stats.put("numPassFilter", 0);
	    stats.put("numFailCoverage", 0);
	    stats.put("numFailFreq", 0);
	    stats.put("numFailQual", 0);
	    stats.put("numFailStrands", 0);
	    stats.put("numFailReads2", 0);
	    stats.put("numFailPvalue", 0);
	    stats.put("numParsingExceptions", 0);
	    stats.put("numNoGenotype", 0);
	    stats.put("numCalledRef", 0);

		// Parse piped input or user-provided pileup file //

	    try
	    {
//	    	 Declare output files //
			PrintStream outFile = null;
			if(params.containsKey("output-file"))
				outFile = new PrintStream( new FileOutputStream(outFileName) );

			PrintStream notFile = null;
			if(params.containsKey("not-file"))
				notFile = new PrintStream( new FileOutputStream(notFileName) );

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
	    			try
	    			{
		    			lineCounter++;

	    				String[] lineContents = line.split("\t");
		    			String chrom = lineContents[0];

		    			if(line.startsWith("#"))
		    			{
		    				// VCF header line //
		    				isVCF = true;
		    				// Print header //
							if(params.containsKey("output-file"))
	    						outFile.println(line);
							else
			    				System.out.println(line);

							if(params.containsKey("not-file"))
								notFile.println(line);
		    			}
		    			else if(chrom.equals("Chrom"))
		    			{
		    				// Native output header line //
							if(params.containsKey("output-file"))
	    						outFile.println(line);
							else
			    				System.out.println(line);

							if(params.containsKey("not-file"))
								notFile.println(line);
		    			}
		    			else
		    			{
		    				int position = Integer.parseInt(lineContents[1]);

		    				// If indel file was provided, check position-1 to position+1 //
		    				boolean indelFilter = false;
		    				if(params.containsKey("indel-file"))
		    				{
		    					String key1 = chrom + "\t" + position;
		    					String key2 = chrom + "\t" + (position - 1);
		    					String key3 = chrom + "\t" + (position + 1);

		    					if(indelPositions.containsKey(key1) || indelPositions.containsKey(key2) || indelPositions.containsKey(key3))
		    					{
		    						indelFilter = true;
		    					}
		    				}

		    				if(indelFilter)
		    				{
		    					stats.put("numNearIndel", (stats.get("numNearIndel") + 1));
			    				stats.put("numVariants", (stats.get("numVariants") + 1));
			    				if(params.containsKey("not-file"))
			    					notFile.println(line);
		    				}
		    				else if(isVCF)
		    				{
		    					int maxCol = lineContents.length;
		    					String vcfLine = "";
		    					int numSamples = 0;
		    					int numSamplesReferencePass = 0;
		    					int numSamplesVariantPass = 0;
		    					for(int colCounter = 0; colCounter < maxCol; colCounter++)
		    					{
		    						if(colCounter < 9)
		    						{
		    							if(colCounter > 0)
		    								vcfLine += "\t";
		    							vcfLine += lineContents[colCounter];
		    						}
		    						else
		    						{
		    							numSamples++;
		    							vcfLine += "\t";
		    							// Evaluate sample //
		    							//GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR//
		    							String[] sampleContents = lineContents[colCounter].split(":");
		    							String gt = sampleContents[0];

		    							stats.put("numVariants", (stats.get("numVariants") + 1));

		    							if(gt.contains("."))
		    							{
		    								// Blank genotype, so ignore //
		    								stats.put("numNoGenotype", (stats.get("numNoGenotype") + 1));
		    							}
		    							else
		    							{
			    							int qualityDepth = Integer.parseInt(sampleContents[3]);
			    							int reads1 = Integer.parseInt(sampleContents[4]);
			    							int reads2 = Integer.parseInt(sampleContents[5]);
			    							double varFreq = (double) reads2 / (double) ((reads1 + reads2));
			    							double pValue = Float.parseFloat(sampleContents[7]);
			    							int qual1 = Integer.parseInt(sampleContents[8]);
			    							int qual2 = Integer.parseInt(sampleContents[9]);
			    							int reads1plus = Integer.parseInt(sampleContents[10]);
			    							int reads1minus = Integer.parseInt(sampleContents[11]);
			    							int reads2plus = Integer.parseInt(sampleContents[12]);
			    							int reads2minus = Integer.parseInt(sampleContents[13]);
			    							boolean strandFail = false;
			    							if(reads1plus > 0 && reads1minus > 0)
			    							{
			    								if(reads2plus == 0 || reads2minus == 0)
			    									strandFail = true;
			    							}

			    							boolean isFiltered = true;

			    							/// Begin checks for either reference or variant //
			    							if(qualityDepth < minCoverage)
			    							{
			    								stats.put("numFailCoverage", (stats.get("numFailCoverage") + 1));
			    								isFiltered = true;
			    							}

			    							if(gt.equals("0/0"))
			    							{
			    								stats.put("numCalledRef", (stats.get("numCalledRef") + 1));

			    								if(!isFiltered)
			    									numSamplesReferencePass++;
			    								// Don't try to filter wild-type calls //
			    							}
			    							else if(reads2 < minReads2)
			    							{
			    								stats.put("numFailReads2", (stats.get("numFailReads2") + 1));
			    								isFiltered = true;
			    							}
			    							else if (qual2 < minAvgQual)
			    							{
			    								stats.put("numFailQual", (stats.get("numFailQual") + 1));
			    								isFiltered = true;
			    							}
			    							else if(varFreq < minVarFreq)
			    							{
			    								stats.put("numFailFreq", (stats.get("numFailFreq") + 1));
			    								isFiltered = true;
			    							}
			    							else if (strandFail)
			    							{
			    								stats.put("numFailStrands", (stats.get("numFailStrands") + 1));
			    								isFiltered = true;
			    							}
			    							else if (pValue > pValueThreshold)
			    							{
			    								stats.put("numFailPvalue", (stats.get("numFailPvalue") + 1));
			    								isFiltered = true;
			    							}
			    							else
			    							{
				    							// Pass the variant //
				    							stats.put("numPassFilter", (stats.get("numPassFilter") + 1));
				    							numSamplesVariantPass++;
			    							}

			    							// IF we have a variant filtered, address that //

			    							if(isFiltered)
			    							{
			    								// Replace genotype with blank value //
			    								lineContents[colCounter].replace(gt, "./.");
			    							}
		    							}


		    							// Append to line //
		    							vcfLine += lineContents[colCounter];
		    						}
		    					}

    							if(numSamplesVariantPass > 0 || numSamplesReferencePass > 0)
    							{
    								if(params.containsKey("output-file"))
    									outFile.println(vcfLine);
    								else
    									System.out.println(vcfLine);
    							}
    							else if(params.containsKey("not-file"))
    							{
    								notFile.println(vcfLine);
    							}

		    				}
		    				else
		    				{
			    				stats.put("numVariants", (stats.get("numVariants") + 1));

			    				// Parse out relevant values //
			    				String ref = lineContents[2];
			    				String var = lineContents[3];
			    				int reads1 = Integer.parseInt(lineContents[4]);
			    				int reads2 = Integer.parseInt(lineContents[5]);
			    				int strands1 = Integer.parseInt(lineContents[7]);
			    				int strands2 = Integer.parseInt(lineContents[8]);
			    				int qual1 = Integer.parseInt(lineContents[9]);
			    				int qual2 = Integer.parseInt(lineContents[10]);
			    				double pValue = Float.parseFloat(lineContents[11]);

			    				int coverage = reads1 + reads2;
			    				double varFreq = (double) reads2 / (double) (reads1 + reads2);

			    				boolean isFiltered = true;


			    				if(coverage >= minCoverage)
			    				{
			    					if(ref.equals(var))
			    					{
			    						stats.put("numCalledRef", (stats.get("numCalledRef") + 1));

    									if(params.containsKey("output-file"))
    			    						outFile.println(line);
    									else
    										System.out.println(line);
			    					}
			    					else if(reads2 >= minReads2)
			    					{
			    						if(strands2 >= minStrands2)
			    						{
			    							if(qual2 >= minAvgQual)
			    							{
			    								if(varFreq >= minVarFreq)
			    								{
			    									// Calculate p-value if it has value of 0.98 meaning not calculated //

			    				    				if(pValue >= 0.98)
			    				    				{
			    				    					pValue = VarScan.getSignificance(reads1, reads2);
			    				    				}


				    								if(pValue <= pValueThreshold)
				    								{
				    									stats.put("numPassFilter", (stats.get("numPassFilter") + 1));

				    									if(params.containsKey("output-file"))
				    			    						outFile.println(line);
				    									else
				    										System.out.println(line);

				    									isFiltered = false;
				    								}
				    								else
				    								{
				    									stats.put("numFailPvalue", (stats.get("numFailPvalue") + 1));
				    								}
			    								}
			    								else
			    								{
			    									stats.put("numFailFreq", (stats.get("numFailFreq") + 1));
			    								}
			    							}
			    							else
		    								{
		    									stats.put("numFailQual", (stats.get("numFailQual") + 1));
		    								}
			    						}
			    						else
	    								{
	    									stats.put("numFailStrands", (stats.get("numFailStrands") + 1));
	    								}
			    					}
			    					else
			    					{
			    						stats.put("numFailReads2", (stats.get("numFailReads2") + 1));
			    					}
			    				}
			    				else
			    				{
			    					stats.put("numFailCoverage", (stats.get("numFailCoverage") + 1));
			    				}

			    				// If not file provided, print this filtered variant to it //

			    				if(isFiltered && params.containsKey("not-file"))
			    					notFile.println(line);


		    				}

		    			}
	    			}
	    			catch(Exception e)
	    			{
	    				if(lineCounter == 1)
	    				{
	    					// Print header //
	    					System.out.println(line);
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

	    		System.err.println(stats.get("numVariants") + " entries in input stream");
	    		System.err.println(stats.get("numNoGenotype") + " had no genotype");
	    		System.err.println(stats.get("numCalledRef") + " were called wild-type");
	    		System.err.println(stats.get("numFailCoverage") + " failed coverage");
	    		System.err.println(stats.get("numFailReads2") + " failed reads2");
	    		System.err.println(stats.get("numFailStrands") + " failed strands");
	    		System.err.println(stats.get("numFailQual") + " failed quality");
	    		System.err.println(stats.get("numFailFreq") + " failed variant frequency < " + minVarFreq);
	    		System.err.println(stats.get("numFailPvalue") + " failed P-value > " + pValueThreshold);
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
    					if(line.startsWith("#"))
    					{
    						// Ignore //
    					}
    					else
    					{
        					String[] lineContents = line.split("\t");
        					String chrom = lineContents[0];
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
}
