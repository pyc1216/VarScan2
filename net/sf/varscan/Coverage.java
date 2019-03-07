/**
 * @(#)FilterSomatic.java
 *
 * Copyright (c) 2009-2010 Daniel C. Koboldt and Washington University in St. Louis
 *
 * COPYRIGHT
 */

package net.sf.varscan;

import java.io.BufferedReader;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.BitSet;
import java.util.HashMap;

/**
 * A class for assessing coverage of target regions (experimental)
 *
 * @version	2.3
 *
 * @author Daniel C. Koboldt <dkoboldt@genome.wustl.edu>
 *
 */
public class Coverage {

	public Coverage(String[] args)
	{
		//		 Define the usage message //
		String usage = "USAGE: java -jar VarScan.jar coverage [pileup-file] OPTIONS\n" +
		"\n" +
		"\tpileup-file - A SAMtools pileup file or piped input\n" +
		"\tOPTIONS:\n" +
		"\t--regions-file\tTab-delimited file of regions of interest (required)\n" +
		"\t--min-base-qual\tMinimum base quality [20]\n" +
		"\t--output-file\tOutput file for coverage report";

		// Set default parameters //

		String outFileName = "";
		String regionsFile = "";
		int minBaseQual = 20;
		int maxDepth = 50;
		String targetFileType = "regions";

		// Parse the input parameters //

		HashMap<String, String> params = VarScan.getParams(args);

		// Adjust parameters based on user input //

		try
		{
			if(params.containsKey("output-file"))
				outFileName = params.get("output-file");

			if(params.containsKey("regions-file"))
				regionsFile = params.get("regions-file");

			if(params.containsKey("min-base-qual"))
				minBaseQual = Integer.parseInt(params.get("min-base-qual"));

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

 	 	// Load target positions //
 	 	HashMap<String, BitSet> targetHash = null;

 	 	if(params.containsKey("regions-file"))
 	 		targetHash = loadTargets(regionsFile, targetFileType, 0);


		// Obtain SAM output with system call to SAMtools //

		try
		{
//			Runtime r = Runtime.getRuntime();
//			Process p = r.exec("samtools view -q " + minMapQual + " " + bamFile);
//			InputStreamReader instream = new InputStreamReader(p.getInputStream());

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

			    	if(numNaps > 10)
			    	{
			    		System.err.println("Input file was not ready after 10 5-second cycles!");
			    		return;
			    	}
		    	}
		    	catch(Exception e)
		    	{

		    	}
	    	}

	    	// Declare an array to count positions at each depth from 0 to max depth //

    	    long[] positionsByDepth = new long[maxDepth + 1];
    	    long basesOnTarget = 0;
    	    long basesOffTarget = 0;

    	    // Prepare to parse the SAM input //

			String line;
			int lineCounter = 0;
			int numParsingExceptions = 0;

	    	if(in != null && in.ready())
	    	{
	    		// Parse the infile line by line //

	    		while ((line = in.readLine()) != null)
	    		{
	    			lineCounter++;//stats.put("numBases", (stats.get("numBases") + 1));

	    			// Output progress line //
	    			if(params.containsKey("verbose") && (lineCounter % 100000) == 0)
		        		System.err.println(lineCounter + " positions parsed...");

	    			// Begin try-catch for line parsing //

	    			try
	    			{
	    				String[] lineContents = line.split("\t");

	    				// Verify expected pileup format //

	    				if(lineContents.length > 5 && lineContents[0].length() > 0 && lineContents[1].length() > 0 && lineContents[2].length() > 0 && lineContents[3].length() > 0)
	    				{
	    					String refName = lineContents[0];
		    	        	int position = Integer.parseInt(lineContents[1]);

		    				// Declare a BitSet //
		    				BitSet refPositions;

		    				boolean inTarget = false;

		    				if(!params.containsKey("regions-file"))
		    				{
		    					// If no regions file provided, report on all positions //
		    					inTarget = true;
		    				}
		    				// Get the position BitSet for this chromosome//
		    				else if(targetHash.containsKey(refName))
		    				{
		    					refPositions = targetHash.get(refName);

		    					// Check to see if position set //
		    					if(refPositions.get(position))
		    					{
		    						inTarget = true;
		    					}
		    				}

		    				if(inTarget)
		    				{
		    					basesOnTarget++;

		    					try {
		    						// Parse out the depth and base qualities //
				    	        	int readDepth = Integer.parseInt(lineContents[3]);
				    	        	String readQualities = lineContents[5];

				    	        	String mapQualities = "";
				    	        	if(lineContents.length > 6)			// Get Map Qualities if available //
				    	        		mapQualities = lineContents[6];

				    	        	int qualDepth = VarScan.qualityDepth(readQualities, minBaseQual);

				    	        	for(int thisDepth = 0; thisDepth <= qualDepth; thisDepth++)
				    	        	{
				    	        		positionsByDepth[thisDepth]++;
				    	        	}
		    					}
		    					catch (Exception e)
		    					{

		    					}

		    				}
		    				else
		    				{
		    					basesOffTarget++;
		    				}


	    				}
	    				else
	    				{
	    					System.err.println("Error: Invalid format for pileup at line " + lineCounter + ":" + line + "\n");
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


		    	// Determine the number of bases targeted //
		    	long targetPositions = 0;
		    	for (String refName : targetHash.keySet())
		    	{
		    		BitSet refPositions = targetHash.get(refName);
		    		for(int thisPosition = 1; thisPosition <= refPositions.size(); thisPosition++)
		    		{
		    			if(refPositions.get(thisPosition))
		    				targetPositions++;
		    		}

		    	}

		    	String outputReport = "";

		    	outputReport += targetPositions + " positions targeted\n";
		    	outputReport += lineCounter + " positions in pileup file\n";

		    	// Print a summary report //
		    	outputReport += "TARGET SPECIFICITY\n";
		    	outputReport += basesOffTarget + " bases off target\n";
		    	outputReport += basesOnTarget + " bases on target\n";

		    	// Print a breadth-by-depth report //

		    	long cov50x = positionsByDepth[50];
		    	long cov40x = positionsByDepth[40] - positionsByDepth[50];
		    	long cov30x = positionsByDepth[30] - positionsByDepth[40];
		    	long cov20x = positionsByDepth[20] - positionsByDepth[30];
		    	long cov10x = positionsByDepth[10] - positionsByDepth[20];
		    	long cov8x = positionsByDepth[8] - positionsByDepth[10];
		    	long cov6x = positionsByDepth[6] - positionsByDepth[8];
		    	long cov2x = positionsByDepth[2] - positionsByDepth[6];
		    	long cov1x = positionsByDepth[1] - positionsByDepth[2];
		    	long cov0x = targetPositions - positionsByDepth[1];

		    	outputReport += "COVERAGE BREADTH-BY-DEPTH\n";
		    	outputReport += "min_depth\t50x\t40x\t30x\t20x\t10x\t8x\t6x\t2x\t1x\t0x\n";
		    	outputReport += "num_positions\t" + cov50x + "\t" + cov40x + "\t" + cov30x + "\t" + cov20x + "\t" + cov10x + "\t" + cov8x + "\t" + cov6x + "\t" + cov2x + "\t" + cov1x + "\t" + cov0x + "\n";

		    	// Print the summary report to error output //
		    	System.err.println(outputReport);

		    	outputReport += "COVERED BASES BY DEPTH\n";
		    	outputReport += "Depth\tBasesCovered\n";
		    	for(int thisDepth = 0; thisDepth <= maxDepth; thisDepth++)
		    	{
		    		outputReport += thisDepth + "\t" + positionsByDepth[thisDepth] + "\n";
		    	}

		    	if(params.containsKey("output-file"))
		    	{
		    		PrintStream outFile = new PrintStream( new FileOutputStream(outFileName) );
		    		outFile.println(outputReport);
		    		outFile.close();
		    	}
	    	}



		}

		catch(Exception e)
		{
	    	System.err.println("Error extracting SAM from BAM: " + e.getLocalizedMessage());
	    	e.printStackTrace(System.err);
	    	return;
		}



	}



	/**
	 * Saves positions into a BitSet hash by chromosome
	 *
	 * @param	fileName	Name of file to be parsed
	 * @param	fileType	Type of file ("positions" or "regions")
	 * @return			HashMap of parameter names and their values
	 */
	static HashMap<String, BitSet> loadTargets(String fileName, String fileType, int marginSize)
	{
		HashMap<String, BitSet> positionsByChrom = new HashMap<String, BitSet>();

		int numRegions = 0;
		int numBases = 0;

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

	    				// Get or create BitSet for this refName //
	    				BitSet refPositions;

	    				if(positionsByChrom.containsKey(refName))
	    				{
	    					refPositions = positionsByChrom.get(refName);
	    				}
	    				else
	    				{
	    					refPositions = new BitSet();
	    				}

	    				// Mark position or regions, depending on what was provided //
	    				int chrStart = 0;
	    				int chrStop = 0;

	    				if(fileType.equals("positions") && lineContents.length > 1)
	    				{
		    				// Set the position to true //
		    				int position = Integer.parseInt(lineContents[1]);
		    				chrStart = position - marginSize;
		    				chrStop = position + marginSize;
	    				}
	    				else if(fileType.equals("regions") && lineContents.length > 2)
	    				{
	    					chrStart = Integer.parseInt(lineContents[1]) - marginSize;
	    					chrStop = Integer.parseInt(lineContents[2]) + marginSize;
	    				}

    					// Check that it won't be an infinite loop//
    					if(chrStart <= chrStop)
    					{
    						numRegions++;

	    					// Mark every position //
	    					for(int position = chrStart; position <= chrStop; position++)
	    					{
	    						if(!refPositions.get(position))
	    						{
	    							numBases++;
	    							refPositions.set(position, true);
	    						}
	    					}
    					}

	    				// Return it to the hash //
	    				positionsByChrom.put(refName, refPositions);
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

		System.err.println(numRegions + " regions parsed");
		System.err.println(numBases + " unique positions targeted");
		return(positionsByChrom);
	}

}
