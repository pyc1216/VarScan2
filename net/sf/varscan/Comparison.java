/**
 * @(#)Comparison.java
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
 * A class for comparing positions (variants) between two files
 *
 * @version	2.3
 *
 * @author Daniel C. Koboldt <dkoboldt@genome.wustl.edu>
 *
 */
public class Comparison {

	public Comparison(String[] args)
	{
		String usage = "USAGE: java -jar VarScan.jar compare [file1] [file2] [type] [output] OPTIONS\n" +
		"\tfile1 - A file of chromosome-positions, tab-delimited\n" +
		"\tfile2 - A file of chromosome-positions, tab-delimited\n" +
		"\ttype - Type of comparison [intersect|merge|unique1|unique2]\n" +
		"\toutput - Output file for the comparison result\n";

		if(args.length < 5)
		{
			 System.out.println(usage);
			 return;
		}

		// Get the required arguments //
		String fileName1 = args[1];
		String fileName2 = args[2];
		String comparisonType = args[3];
		String outFileName = args[4];

		System.err.println("File 1: " + fileName1);
		System.err.println("File 2: " + fileName2);

		try
		{
			// Declare output file //
			PrintStream outFile = null;

			outFile = new PrintStream( new FileOutputStream(outFileName) );

	 		BufferedReader file1 = new BufferedReader(new FileReader(fileName1));
	 		BufferedReader file2 = new BufferedReader(new FileReader(fileName2));

		    if(!(file1.ready() && file2.ready()))
		    {
		    	System.err.println("ERROR: Invalid input file(s)");
		    	return;
		    }

		    // Load the positions of both files into a hash //
		    System.err.println("Loading positions from file 1");
		    HashMap<String, BitSet> positionHash1 = loadPositions(fileName1);
		    System.err.println("Loading positions from file 2");
		    HashMap<String, BitSet> positionHash2 = loadPositions(fileName2);
		    System.err.println("Done");

		    // Reset counters ///

		    int numShared = 0;
		    int uniqueToFile1 = 0;
		    int uniqueToFile2 = 0;

		    // Parse the lines in file 1 //
		    String line = "";
		    int lineCounter = 0;

		    while ((line = file1.readLine()) != null)
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

	    				// Declare booleans //

	    				boolean inFile1 = false;
	    				boolean inFile2 = false;

	    				// Declare a BitSet //
	    				BitSet refPositions;

	    				if(positionHash1.containsKey(refName))
	    				{
	    					refPositions = positionHash1.get(refName);
	    					if(refPositions.get(position))
	    						inFile1 = true;
	    				}

	    				if(positionHash2.containsKey(refName))
	    				{
	    					refPositions = positionHash2.get(refName);
	    					if(refPositions.get(position))
	    						inFile2 = true;
	    				}

	    				// Check to see if shared //
	    				if(inFile1 && inFile2)
	    				{
	    					numShared++;
	    					if(comparisonType.equals("intersect"))
	    					{
	    						outFile.println(line);
	    					}
	    				}
	    				else if(inFile1)
	    				{
	    					if(comparisonType.equals("unique1"))
	    						outFile.println(line);
	    					uniqueToFile1++;
	    				}

	    				// Check to see if merging //
	    				if(comparisonType.equals("merge"))
	    				{
	    					outFile.println(line);
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


		    while ((line = file2.readLine()) != null)
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

	    				// Declare booleans //

	    				boolean inFile1 = false;
	    				boolean inFile2 = false;

	    				// Declare a BitSet //
	    				BitSet refPositions;

	    				if(positionHash1.containsKey(refName))
	    				{
	    					refPositions = positionHash1.get(refName);
	    					if(refPositions.get(position))
	    						inFile1 = true;
	    				}

	    				if(positionHash2.containsKey(refName))
	    				{
	    					refPositions = positionHash2.get(refName);
	    					if(refPositions.get(position))
	    						inFile2 = true;
	    				}

	    				// Check to see if shared //
	    				if(inFile1 && inFile2)
	    				{
	    					// Already counted and printed in file 1 //
	    				}
	    				else if(inFile2)
	    				{
	    					if(comparisonType.equals("merge") || comparisonType.equals("unique2"))
	    						outFile.println(line);

	    					uniqueToFile2++;
	    				}

	    			}
	    			catch(Exception e)
	    			{
	    				if(lineCounter == 1)
	    					outFile.println(line);
	    				else
	    					System.err.println("Warning: Unable to parse chrom/position from " + line);
	    			}
	    		}
		    }

		    file1.close();
		    file2.close();

		    int numTotal = numShared + uniqueToFile1 + uniqueToFile2;
		    System.err.println(numTotal + " total positions");
		    System.err.println(uniqueToFile1 + " positions unique to file 1");
		    System.err.println(uniqueToFile2 + " positions unique to file 2");
		    System.err.println(numShared + " positions shared");
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
	 * @param	args	Command-line arguments
	 * @return			HashMap of parameter names and their values
	 */
	static HashMap<String, BitSet> loadPositions(String fileName)
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
	    		if(lineContents.length >= 2)
	    		{
	    			// Try to parse chrom and position //
	    			try
	    			{
	    				String refName = lineContents[0];
	    				int position = Integer.parseInt(lineContents[1]);

	    				// Get or create BitSet //
	    				BitSet refPositions;

	    				if(positionsByChrom.containsKey(refName))
	    				{
	    					refPositions = positionsByChrom.get(refName);
	    				}
	    				else
	    				{
	    					refPositions = new BitSet(position + 1);
	    				}

	    				// Set the position to true //
	    				refPositions.set(position, true);

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
