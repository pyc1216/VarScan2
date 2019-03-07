/**
 * @(#)ReadCounts.java
 *
 * Copyright (c) 2009-2010 Daniel C. Koboldt and Washington University in St. Louis
 *
 * COPYRIGHT
 */

package net.sf.varscan;

import java.io.BufferedReader;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.BitSet;
import java.util.HashMap;

/**
 * A class for obtaining read counts for a list of variants
 *
 * @version	2.3
 *
 * @author Daniel C. Koboldt <dkoboldt@genome.wustl.edu>
 *
 */
public class ReadCounts {

	public ReadCounts(String[] args, HashMap<String, String> params)
	{
		String usage = "USAGE: java -jar VarScan.jar readcounts [pileup] OPTIONS\n" +
		"\tOPTIONS:\n" +
		"\t--variants-file\tA list of variants at which to report readcounts\n" +
		"\t--output-file\tOutput file to contain the readcounts\n" +
		"\t--min-coverage\tMinimum read depth at a position to make a call [1]\n" +
		"\t--min-base-qual\tMinimum base quality at a position to count a read [20]\n";

		// Set parameter defaults //

		String variantsFile = "";
		String outputFile = "";
		int minCoverage = 1;
		int minBaseQual = 20;

		// Get any user-provided parameters //

		try
		{
			if(params.containsKey("min-coverage"))
				 minCoverage = Integer.parseInt(params.get("min-coverage"));

			if(params.containsKey("min-base-qual"))
				 minBaseQual = Integer.parseInt(params.get("min-base-qual"));

			if(params.containsKey("variants-file"))
				 variantsFile = params.get("variants-file");

			if(params.containsKey("output-file"))
				 outputFile = params.get("output-file");

			 System.err.println("Min coverage:\t" + minCoverage);
			 System.err.println("Min base qual:\t" + minBaseQual);
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

	    // Define the statistics hash and reset counters //

	    HashMap<String, Integer> stats = new HashMap<String, Integer>();
	    stats.put("numPositions", 0);
	    stats.put("numIncluded", 0);
	    stats.put("numCovered", 0);

		//	If a list of variants was provided, save it //

		HashMap<String, String> variantPositions = new HashMap<String, String>();
		if(params.containsKey("variants-file"))
		{
			System.err.println("Loading variant positions from " + variantsFile);
			variantPositions = loadVariants(variantsFile);
			System.err.println(variantPositions.size() + " variant positions saved");
		}

		// Attempt to obtain and parse pileup input //

	    try
	    {
	    	// If output file was provided, open it //
	    	PrintStream out = null; // declare a print stream object

	    	if(params.containsKey("output-file"))
	    	{
	    		out = new PrintStream( new FileOutputStream(outputFile) );
	    		out.println("chrom\tposition\tref_base\tdepth\tq" + minBaseQual + "_depth\tbase:reads:strands:avg_qual:map_qual:plus_reads:minus_reads");
	    	}

	    	// Declare file-parsing variables //

	    	BufferedReader in = VarScan.getInfile(args);
	    	String line;

	    	// If no input, print usage //

	    	if(in == null)
	    	{
	    		System.out.println(usage);
	    		System.exit(10);
	    	}

	    	// If input file not ready, give it a few seconds //
	    	int numNaps = 0;

	    	if(!in.ready())
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

		    	}
	    	}

	    	// If pileup input was provided, begin parsing it //

	    	if(in.ready())
	    	{
	    		while ((line = in.readLine()) != null)
	    		{
	    			stats.put("numPositions", (stats.get("numPositions") + 1));

	    			// Output progress line //
	    			if(params.containsKey("verbose") && (stats.get("numPositions") % 100000) == 0)
		        		System.err.println(stats.get("numPositions") + " positions parsed...");

	    			// Begin try-catch for line parsing //

	    			try
	    			{
	    				String[] lineContents = line.split("\t");

	    				// Verify expected pileup format //

	    				if(lineContents.length > 5 && lineContents[0].length() > 0 && lineContents[1].length() > 0 && lineContents[2].length() > 0 && lineContents[3].length() > 0)
	    				{
	    					String refName = lineContents[0];
		    	        	String position = lineContents[1];
		    	        	String refBase = lineContents[2].toUpperCase();
		    	        	int readDepth = Integer.parseInt(lineContents[3]);
		    	        	String readBases = lineContents[4];
		    	        	String readQualities = lineContents[5];
		    	        	String mapQualities = "";
		    	        	if(lineContents.length > 6)			// Get Map Qualities if available //
		    	        		mapQualities = lineContents[6];

		    	        	// If variant file was provided, verify that this position matches one in list //

		    	        	if(!params.containsKey("variants-file") || variantPositions.containsKey(refName + "\t" + position))
		    	        	{
		    	        		stats.put("numIncluded", (stats.get("numIncluded") + 1));

		    	        		// Build output line //
		    	        		String outputLine = refName + "\t" + position + "\t" + refBase + "\t" + readDepth + "\t";

			    	        	if(readDepth >= minCoverage)
			    	        	{
			    	        		stats.put("numCovered", (stats.get("numCovered") + 1));

			    	        		// Obtain the readcounts //

			    	        		HashMap<String, String> readCounts = VarScan.getReadCounts(refBase, readBases, readQualities, minBaseQual, mapQualities);

			    	        		// Build array of allele keys observed and sort it //
			    	        		String[] alleleKeys = (String[]) readCounts.keySet().toArray(new String[0]);
			    	    			Arrays.sort(alleleKeys);

			    	    			// Calculate the # of reads that met quality threshold //
			    	    			int readDepthQual = VarScan.qualityDepth(readQualities, minBaseQual);


			    	    			// Add the quality-met read count to output //

			    	    			outputLine += readDepthQual + "\t";

			    	        		// First, get the ref result //
			    	    			//reads2 strands2 qual2 map2 readsPlus readsMinus
			    	        		String refResult = "0\t0\t0\t0\t0\t0";
			    	        		if(readCounts.containsKey(refBase))
			    	        			refResult = readCounts.get(refBase);

			    	        		// Replace tabs with colons for ref base //
			    	        		refResult = refBase + ":" + refResult.replace("\t", ":");

			    	        		outputLine += refResult + "\t";

			    	        		// If a variant file was provided, try to get the desired variant allele //

			    	        		String desiredAllele = "";
			    	        		if(variantPositions.containsKey(refName + "\t" + position))
			    	        		{
			    	        			try
			    	        			{
				    	        			String[] varContents = variantPositions.get(refName + "\t" + position).split("\t");
				    	        			if(varContents.length > 1 && varContents[1].length() > 0)
				    	        			{
				    	        				desiredAllele = varContents[1];

				    	        				String varResult = "0\t0\t0\t0\t0\t0";
						    	        		if(readCounts.containsKey(desiredAllele))
						    	        			varResult = readCounts.get(desiredAllele);

						    	        		outputLine += desiredAllele + "\t" + varResult + "\t";
				    	        			}
			    	        			}
			    	        			catch(Exception e)
			    	        			{
			    	        				System.err.println("Warning: Error parsing variant position entry: " + variantPositions.get(refName + "\t" + position) + "\n" + e.getLocalizedMessage());
			    	        		    	e.printStackTrace(System.err);
			    	        			}
			    	        		}

			    	        		// Go through all bases observed //

			    	    			for(String allele : alleleKeys)
			    	    			{
			    	    				String[] alleleContents = readCounts.get(allele).split("\t");
			    	    				if(allele.equals(refBase) || allele.equals(desiredAllele))
			    	    				{
			    	    					// Skip the reference base and desired base //
			    	    				}
			    	    				else
			    	    				{
			    	    					try {
				    	    					int thisReads2 = Integer.parseInt(alleleContents[0]);
				    	    					int thisStrands2 = Integer.parseInt(alleleContents[1]);
				    	    					int thisAvgQual2 = Integer.parseInt(alleleContents[2]);
				    	    					int thisMapQual2 = Integer.parseInt(alleleContents[3]);
				    	    					int thisReads2plus = Integer.parseInt(alleleContents[4]);
				    	    					int thisReads2minus = Integer.parseInt(alleleContents[5]);

				    	    					String varResult = allele + ":" + thisReads2 + ":" + thisStrands2 + ":" + thisAvgQual2 + ":" + thisMapQual2 + ":" + thisReads2plus + ":" + thisReads2minus;
				    	    					outputLine += varResult + "\t";
			    	    					}
			    	    					catch(Exception e)
			    	    					{

			    	    					}

			    	    				}
			    	    			}

			    	    			// Print the output line //

			    	    			if(params.containsKey("output-file"))
			    	    			{
			    	    				out.println(outputLine);
			    	    			}
			    	    			else
			    	    			{
			    	    				System.err.println(outputLine);
			    	    			}
			    	        	}
			    	        	else if(variantPositions.containsKey(refName + "\t" + position))
			    	        	{
			    	    			// Insufficient coverage - Print the output line //

			    	    			if(params.containsKey("output-file"))
			    	    			{
			    	    				out.println(outputLine);
			    	    			}
			    	    			else
			    	    			{
			    	    				System.err.println(outputLine);
			    	    			}
			    	        	}
		    	        	}

	    				}
	    				else
	    				{
	    					System.err.println("Error: Invalid format for pileup at line " + stats.get("numBases") + "\n" + line + "\n");
	    					return;
	    				}
	    			}
	    			catch(Exception e)
	    		    {
	    		    	System.err.println("Parsing Exception on line:\n" + line + "\n" + e.getLocalizedMessage());
	    				stats.put("numParsingExceptions", (stats.get("numParsingExceptions") + 1));
	    				if(stats.get("numParsingExceptions") >= 5)
	    				{
	    					System.err.println("Too many parsing exceptions encountered; exiting");
	    					return;
	    				}
	    		    	return;
	    		    }


	    		}
	    	}
	    	else
	    	{
	    		System.err.println("Input was not ready for parsing!");
		    	return;
	    	}

	    	in.close();

	    	if(params.containsKey("output-file"))
	    	{
	    		out.close();
	    	}

	    	// Print summary statistics //

	    	System.err.println(stats.get("numPositions") + " positions in pileup file");
	    	System.err.println(stats.get("numIncluded") + " included in readcount analysis");
	    	System.err.println(stats.get("numCovered") + " met minimum coverage");
	    }
	    catch(Exception e)
	    {
	    	System.err.println("Error parsing input: " + e.getLocalizedMessage());
	    	e.printStackTrace(System.err);
	    	System.exit(11);
	    }

	}


	/**
	 * Saves variants into a HashMap
	 *
	 * @param	args	Command-line arguments
	 * @return			HashMap of parameter names and their values
	 */
	static HashMap<String, String> loadVariants(String fileName)
	{
		HashMap<String, String> variants = new HashMap<String, String>();

		try
		{
			BufferedReader infile = new BufferedReader(new FileReader(fileName));

			String line = "";
			int lineCounter = 0;

	    	while ((line = infile.readLine()) != null)
	    	{
	    		lineCounter++;

	    		String[] lineContents = line.split("\t");
	    		if(lineContents.length >= 2)
	    		{
	    			// Try to parse chrom and position //
	    			try
	    			{
	    				String refName = lineContents[0];
	    				int position = Integer.parseInt(lineContents[1]);

	    				String allele1 = "";
	    				String allele2 = "";

	    				try{
	    					allele1 = lineContents[2];
	    					allele2 = lineContents[3];
	    				}
	    				catch(Exception e)
	    				{
	    					// Don't load variant positions //
	    				}

	    				String positionKey = refName + "\t" + position;
	    				variants.put(positionKey, allele1 + "\t" + allele2);
	    			}
	    			catch(Exception e)
	    			{
	    				if(lineCounter > 1)
	    					System.err.println("Warning: Unable to parse chrom/position from " + line);
	    			}


	    		}
	    	}

	    	infile.close();
		}
		catch(Exception e)
		{
			System.err.println("ERROR: File Parsing Exception: " + e.getLocalizedMessage());
			e.printStackTrace(System.err);
		}


		return(variants);
	}
}
