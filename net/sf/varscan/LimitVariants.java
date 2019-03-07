/**
 * @(#)LimitVariants.java
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
import java.util.HashMap;
import java.util.BitSet;

/**
 * A class for restricting a list of variants to a given set of positions or regions
 *
 * @version	2.3
 *
 * @author Daniel C. Koboldt <dkoboldt@genome.wustl.edu>
 *
 */
public class LimitVariants {

	public LimitVariants(String[] args)
	{
		String usage = "USAGE: java -jar VarScan.jar limit [infile] OPTIONS\n" +
		"\tinfile - A file of chromosome-positions, tab-delimited\n" +
		"\tOPTIONS\n" +
		"\t--positions-file - a file of chromosome-positions, tab delimited, or VCF\n" +
		"\t--regions-file - a file of chromosome-start-stops, tab delimited\n" +
		"\t--margin-size - shoulder bases to allow on either side of targets [0]\n" +
		"\t--output-file - Output file for the matching variants\n" +
		"\t--not-file - Output file for variants NOT matching regions/positions\n";

		// Declare argument variables //
		String outFileName = "";
		String notFileName = "";
		String targetFileName = "";
		String targetFileType = "";
		int marginSize = 0;

		// Parse command-line parameters //
		HashMap<String, String> params = VarScan.getParams(args);

		// Try adjusting any provided parameters based on user inut //

		try
		{
			if(params.containsKey("output-file"))
				outFileName = params.get("output-file");

			if(params.containsKey("not-file"))
				notFileName = params.get("not-file");


			if(params.containsKey("positions-file"))
			{
				targetFileName = params.get("positions-file");
				targetFileType = "positions";
			}
			else if(params.containsKey("regions-file"))
			{
				targetFileName = params.get("regions-file");
				targetFileType = "regions";
			}
			else
			{
				System.err.println("Please provide a regions file or a positions file");
				System.err.println(usage);
				return;
			}

			if(params.containsKey("margin-size"))
			{
				marginSize = Integer.parseInt(params.get("margin-size"));
			}
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

		try
		{
	 	 	// Obtain input file //
			BufferedReader infile = VarScan.getInfile(args);

			// Declare output file //
			PrintStream outFile = null;
			if(params.containsKey("output-file"))
				outFile = new PrintStream( new FileOutputStream(outFileName) );

			// Declare not output file //
			PrintStream notFile = null;
			if(params.containsKey("not-file"))
				notFile = new PrintStream( new FileOutputStream(notFileName) );

	 	 	// Load target positions //
	 	 	HashMap<String, BitSet> targetHash = loadTargets(targetFileName, targetFileType, marginSize);

	 	 	// Declare file-parsing variables //
	 	 	String line = "";
	 	 	int lineCounter = 0;
	 	 	int numVariants = 0;
	 	 	int numInTarget = 0;

	 	 	// Actually parse the infile //

	 	 	while ((line = infile.readLine()) != null)
		    {
		    	lineCounter++;

	    		String[] lineContents = line.split("\t");
	    		if(line.substring(0, 1).equals("#"))
	    		{
	    			// Handle VCF headers //
	    			outFile.println(line);
	    			notFile.println(line);
	    		}
	    		else if(line.substring(0, 5).toLowerCase().equals("chrom"))
	    		{
	    			// Handle native file headers //
	    			outFile.println(line);
	    			if(params.containsKey("not-file"))
	    				notFile.println(line);
	    		}
	    		else if(lineContents.length >= 2)
	    		{
	    			// Try to parse chrom and position //
	    			try
	    			{
	    				String refName = lineContents[0];
	    				int position = Integer.parseInt(lineContents[1]);

	    				numVariants++;

	    				// Declare a BitSet //
	    				BitSet refPositions;

	    				boolean inTarget = false;

	    				// Get the position BitSet for this chromosome//
	    				if(targetHash.containsKey(refName))
	    				{
	    					refPositions = targetHash.get(refName);

	    					// Check to see if position set //
	    					if(refPositions.get(position))
	    					{
	    						inTarget = true;
    							numInTarget++;
	    						if(params.containsKey("output-file"))
	    						{
	    							outFile.println(line);
	    						}
	    					}
	    				}

	    				// If no match and not file declared, print to it //
	    				if(!inTarget && params.containsKey("not-file"))
	    				{
	    					notFile.println(line);
	    				}

	    			}
	    			catch(Exception e)
	    			{
	    				if(lineCounter == 1)
	    				{
//	    					 Skip header // outFile.println(line);
	    				}
	    				else
	    					System.err.println("Warning: Unable to parse chrom/position from " + line);
	    			}
	    		}
		    }

	 	 	float pctInTarget = (float) 0;

	 	 	if(numVariants > 0 && numInTarget > 0)
	 	 	{
	 	 		pctInTarget = (float) numInTarget / (float) numVariants * (float) 100;
	 	 	}

	 	 	// Print summary statistics //

	 	 	System.err.println(numVariants + " variants in input file");
	 	 	System.err.println(numInTarget + " variants (" + pctInTarget + "%) matched target positions");
		}
		catch(Exception e)
		{
			System.err.println("ERROR: File Parsing Exception: " + e.getLocalizedMessage());
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

		try
		{
			BufferedReader infile = new BufferedReader(new FileReader(fileName));

			String line = "";
			int lineCounter = 0;

	    	while ((line = infile.readLine()) != null)
	    	{
	    		lineCounter++;
	    		String[] lineContents = line.split("\t");

	    		if(line.substring(0, 1).equals("#"))
	    		{
	    			// Ignore VCF headers //
	    		}
	    		else if(lineContents.length >= 2)
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
	    					// Mark every position //
	    					for(int position = chrStart; position <= chrStop; position++)
	    					{
	    						refPositions.set(position, true);
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


		return(positionsByChrom);
	}

}
