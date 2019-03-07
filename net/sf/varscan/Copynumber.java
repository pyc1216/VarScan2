package net.sf.varscan;

import java.io.BufferedReader;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.HashMap;
/**
 * A class for calling copy number variants between a tumor and a matched normal sample
 *
 * @version	2.3
 *
 * @author Daniel C. Koboldt <dkoboldt@genome.wustl.edu>
 *
 */
public class Copynumber {
	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Constructor with two arguments (string[], boolean) expects mpileup input 					  //
	////////////////////////////////////////////////////////////////////////////////////////////////////
	public Copynumber(String[] args, boolean isMpileup)
	{
		String usage = "USAGE: java -jar VarScan.jar copynumber [normal-tumor.mpileup] [Opt: output] OPTIONS\n" +
			"\tnormal-tumor.mpileup - The SAMtools mpileup file for Normal and Tumor\n" +
			"\toutput - Output base name for files\n" +
			"\nOPTIONS:\n" +
			"\t--min-base-qual - Minimum base quality to count for coverage [20]\n" +
			"\t--min-map-qual - Minimum read mapping quality to count for coverage [20]\n" +
			"\t--min-coverage - Minimum coverage threshold for copynumber segments [20]\n" +
			"\t--min-segment-size - Minimum number of consecutive bases to report a segment [10]\n" +
			"\t--max-segment-size - Max size before a new segment is made [100]\n" +
			"\t--p-value - P-value threshold for significant copynumber change-point [0.01]\n" +
			"\t--data-ratio - The normal/tumor input data ratio for copynumber adjustment [1.0]\n";

		if(args.length < 2)
		{
			 System.err.println(usage);
			 return;
		}

		// Set parameter defaults //

		HashMap<String, String> params = VarScan.getParams(args);

		// Set up formatting for p-values //
		DecimalFormat pvalueFormat = new DecimalFormat("0.####E0");

		String outputName = "output";

		if(args.length >= 3 && !args[2].startsWith("-"))
		{
			outputName = args[2];
		}

		//	Set parameter defaults //

		int minCoverage = 	10;
		int minBaseQual = 	15;
		int minSegmentSize = 10;
		int maxSegmentSize = 100;
		double dataRatio = 1.00;
		double pValueThreshold = 0.01;

		// Try adjusting any provided parameters based on user inut //
		try
		{
			if(params.containsKey("min-coverage"))
			{
				 minCoverage = Integer.parseInt(params.get("min-coverage"));
			}

			if(params.containsKey("min-base-qual"))
				 minBaseQual = Integer.parseInt(params.get("min-base-qual"));

			if(params.containsKey("min-segment-size"))
				 minSegmentSize = Integer.parseInt(params.get("min-segment-size"));

			if(params.containsKey("max-segment-size"))
				 maxSegmentSize = Integer.parseInt(params.get("max-segment-size"));

			if(params.containsKey("p-value"))
				 pValueThreshold = Double.parseDouble(params.get("p-value"));

			if(params.containsKey("data-ratio"))
				 dataRatio = Double.parseDouble(params.get("data-ratio"));

			System.err.println("Min coverage:\t" + minCoverage);
			System.err.println("Min avg qual:\t" + minBaseQual);
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

		// Check for correct input //

		if(args.length < 3)
		{
			System.err.println("Please provide an output file basename!");
			System.err.println(usage);
			System.exit(1);
		}


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
				// Declare output file //
		 	 	PrintStream outCopySegments = null; // declare a print stream object for copynumber segments

		 		outCopySegments = new PrintStream( new FileOutputStream(outputName + ".copynumber") );
		 		outCopySegments.println("chrom\tchr_start\tchr_stop\tnum_positions\tnormal_depth\ttumor_depth\tlog2_ratio\tgc_content");


	    		System.err.println("Reading mpileup input...");
	    		int numParsingExceptions = 0;

	    		// Statistics counters //
	    		long sharedPositions = 0;
	    		long comparedPositions = 0;
	    		long rawCopySegments = 0;
	    		long goodCopySegments = 0;

	    		// Set some default parsing variables //
			    String chromNormal = "";
			    String chromTumor = "";
			    String refBase = "";
			    int posNormal = 0;
			    int posTumor = 0;

			    // Parameters for copy number calling //
			    String copyChrom = "";
			    int copyStart = 0;
			    int copyStop = 0;
			    int copyDepthNormal = 0;
			    int copyDepthTumor = 0;
			    long copySumNormal = 0;
			    long copySumTumor = 0;
			    long copyPositions = 0;
			    long copyPositionsGC = 0;

			    DecimalFormat oneDigit = new DecimalFormat("#0.0");
			    DecimalFormat threeDigits = new DecimalFormat("#0.000");

	    		// Parse the infile line by line //

	    		while ((line = in.readLine()) != null)
	    		{

	    			// Begin try-catch for line parsing //

	    			try
	    			{
	    				String[] lineContents = line.split("\t");

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
		    	        	//String normalBases = lineContents[normalOffset + 1];
		    	        	String normalQualities = lineContents[normalOffset + 2];

	    					// Parse tumor, which should be second sample //
	    					int tumorOffset = 6;
	    					int pileupDepthTumor = Integer.parseInt(lineContents[tumorOffset]);
		    	        	//String tumorBases = lineContents[tumorOffset + 1];
		    	        	String tumorQualities = lineContents[tumorOffset + 2];


		    	        	// If either sample met the minimum coverage and both had at least one read //

//		    	        	if((pileupDepthNormal >= minCoverage || pileupDepthTumor >= minCoverage) && normalQualities.length() > 0)// && tumorQualities.length() > 0)

		    	        	// We want the normal sample to meet the minimum coverage because that's the comparator //
		    	        	if(pileupDepthNormal >= minCoverage && normalQualities.length() > 0)// && tumorQualities.length() > 0)
	    					{
	    						comparedPositions++;
	    						// Get the depth of bases above minimum quality //

	    	    				int normalDepth = VarScan.qualityDepth(normalQualities, minBaseQual);
	    	    				int tumorDepth = 0;
	    	    				if(tumorQualities.length() > 0)
	    	    					tumorDepth = VarScan.qualityDepth(tumorQualities, minBaseQual);

	    	    				// Determine if we have a copy changepoint //
	    	    				// If this base is not contiguous with the copyRegion
	    	    				// If the normal or tumor depth changes //

	    	    				int diffNormal = Math.abs(copyDepthNormal - normalDepth);
	    	    				int diffTumor = Math.abs(copyDepthTumor - tumorDepth);
	    	    				int posDiff = posTumor - copyStop;

	    	    				// DETERMINE IF WE CONTINUE THIS REGION OR PROCESS IT AND START A NEW ONE //

	    	    				boolean continueFlag = false;

	    	    				// If chromosomes differ or contiguity broken, process the region //

	    	    				if(posDiff > 2 || !(copyChrom.equals(chromTumor)))
	    	    				{
	    	    					continueFlag = false;
	    	    				}
	    	    				else
	    	    				{
	    	    					if(copyPositions >= maxSegmentSize)
	    	    					{
	    	    						continueFlag = false;
	    	    					}
	    	    					else if(diffNormal <= 2 && diffTumor <= 2)
	    	    					{
	    	    						continueFlag = true;
	    	    					}
	    	    					else
	    	    					{
	    	    						// Do a Fisher's exact test on the copy number changes. ##

	        	    					double changePvalue = VarScan.getSignificance(copyDepthNormal, copyDepthTumor, normalDepth, tumorDepth);

	        	    					// If depth change not significant, continue with region //
	        	    					if(changePvalue >= pValueThreshold)
	        	    					{
	        	    						continueFlag = true;
	        	    					}
	        	    					else
	        	    					{
	        	    						continueFlag = false;
	        	    					}

	    	    					}
	    	    				}


	    	    				// If continuing, extend this region and don't process yet //

	    	    				if(continueFlag)
	    	    				{
	    	    					copySumNormal += normalDepth;
	    	    					copySumTumor += tumorDepth;
	    	    					copyPositions++;
	    	    					if(refBase.equals("G") || refBase.equals("C") || refBase.equals("g") || refBase.equals("c"))
	    	    						copyPositionsGC++;
	    	    					copyStop = posTumor;
	    	    				}

	    	    				// Otherwise, process this region (if it qualifies) and start a new one //

	    	    				else
	    	    				{
	    	    					if(copyPositions >= minSegmentSize)
	    	    					{
	    	    						rawCopySegments++;
	    	    						String regionResults = processCopyRegion(copyChrom, copyStart, copyStop, copyPositions, copyPositionsGC, copySumNormal, copySumTumor, minCoverage, dataRatio);

	    	    						if(regionResults.length() > 0)
	    	    						{
	    	    							outCopySegments.println(regionResults);
	    	    							goodCopySegments++;
	    	    						}
	    	    					}

	    	    					// Start a new copyNumber region //
	    	    					copyChrom = chromTumor;
	    	    					copyStart = posTumor;
	    	    					copyStop = posTumor;
	    	    					copyDepthNormal = normalDepth;
	    	    					copyDepthTumor = tumorDepth;
	    	    					copySumNormal = normalDepth;
	    	    					copySumTumor = tumorDepth;
	    	    					copyPositions = 1;
	    	    					if(refBase.equals("G") || refBase.equals("C") || refBase.equals("g") || refBase.equals("c"))
	    	    						copyPositionsGC = 1;
	    	    					else
	    	    						copyPositionsGC = 0;
	    	    				}


	    					}
	    					else
	    					{
	    						// If minimum coverage was not met, print region //
		    					// If we had a copyNumber region that met minimum coverage, report it //
		    					if(copyPositions >= minSegmentSize)
		    					{
		    						rawCopySegments++;
		    						String regionResults = processCopyRegion(copyChrom, copyStart, copyStop, copyPositions, copyPositionsGC, copySumNormal, copySumTumor, minCoverage, dataRatio);

		    						if(regionResults.length() > 0)
		    						{
		    							outCopySegments.println(regionResults);
		    							goodCopySegments++;
		    						}
		    					}

		    					// Reset the copyNumber region //
		    					copyChrom = "";
		    					copyStart = 0;
		    					copyStop = 0;
		    					copyDepthNormal = 0;
		    					copyDepthTumor = 0;
		    					copySumNormal = 0;
		    					copySumTumor = 0;
		    					copyPositions = 0;
		    					copyPositionsGC = 0;
	    					}

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

				// Last region: If minimum coverage was not met, print region //
				// If we had a copyNumber region that met minimum coverage, report it //
				if(copyPositions > minSegmentSize)
				{
					rawCopySegments++;
					String regionResults = processCopyRegion(copyChrom, copyStart, copyStop, copyPositions, copyPositionsGC, copySumNormal, copySumTumor, minCoverage, dataRatio);

					if(regionResults.length() > 0)
					{
						outCopySegments.println(regionResults);
						goodCopySegments++;
					}
				}

				in.close();

			    System.err.println(sharedPositions + " positions in mpileup"); //stats.get("sharedPositions")
			    System.err.println(comparedPositions + " had sufficient coverage for comparison"); //stats.get("comparedPositions")
			    System.err.println(rawCopySegments + " raw copynumber segments with size > " + minSegmentSize);
			    System.err.println(goodCopySegments + " good copynumber segments with depth > " + minCoverage);

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
	// Constructor with one argument (string[]) expects independent normal and tumor pileups as input //
	////////////////////////////////////////////////////////////////////////////////////////////////////

	public Copynumber(String[] args)
	{
		String usage = "USAGE: VarScan copynumber [normal_pileup] [tumor_pileup] [Opt: output] OPTIONS\n" +
			"\tnormal_pileup - The SAMtools pileup file for Normal\n" +
			"\ttumor_pileup - The SAMtools pileup file for Tumor\n" +
			"\toutput - Output base name for files\n" +
			"***If you have a single mpileup, see VarScan copynumber -mpileup 1 -h ***\n" +
			"\nOPTIONS:\n" +
			"\t--min-base-qual - Minimum base quality to count for coverage [20]\n" +
			"\t--min-map-qual - Minimum read mapping quality to count for coverage [20]\n" +
			"\t--min-coverage - Minimum coverage threshold for copynumber segments [20]\n" +
			"\t--min-segment-size - Minimum number of consecutive bases to report a segment [10]\n" +
			"\t--max-segment-size - Max size before a new segment is made [100]\n" +
			"\t--p-value - P-value threshold for significant copynumber change-point [0.01]\n" +
			"\t--data-ratio - The normal/tumor input data ratio for copynumber adjustment [1.0]\n";

		if(args.length < 3)
		{
			 System.err.println(usage);
			 return;
		}

		// Get the required arguments //
		String normalPileupFile = args[1];
		String tumorPileupFile = args[2];

		String outputName = "output";

		if(args.length >= 4 && !args[3].startsWith("-"))
		{
			outputName = args[3];
		}

		System.err.println("Normal Pileup: " + normalPileupFile);
		System.err.println("Tumor Pileup: " + tumorPileupFile);

		//	Set parameter defaults //

		int minCoverage = 	10;
		int minBaseQual = 	15;
		int minSegmentSize = 10;
		int maxSegmentSize = 100;
		double dataRatio = 1.00;
		double pValueThreshold = 0.01;

		// Parse command-line parameters //
		HashMap<String, String> params = VarScan.getParams(args);

		// Try adjusting any provided parameters based on user inut //
		try
		{
			if(params.containsKey("min-coverage"))
			{
				 minCoverage = Integer.parseInt(params.get("min-coverage"));
			}

			if(params.containsKey("min-base-qual"))
				 minBaseQual = Integer.parseInt(params.get("min-base-qual"));

			if(params.containsKey("min-segment-size"))
				 minSegmentSize = Integer.parseInt(params.get("min-segment-size"));

			if(params.containsKey("max-segment-size"))
				 maxSegmentSize = Integer.parseInt(params.get("max-segment-size"));

			if(params.containsKey("p-value"))
				 pValueThreshold = Double.parseDouble(params.get("p-value"));

			if(params.containsKey("data-ratio"))
				 dataRatio = Double.parseDouble(params.get("data-ratio"));

			System.err.println("Min coverage:\t" + minCoverage);
			System.err.println("Min avg qual:\t" + minBaseQual);
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

		// Check for correct input //

		if(args.length < 3)
		{
			System.err.println("Please provide an output file basename!");
			System.err.println(usage);
			System.exit(1);
		}


		// Statistics counters //
		long tumorPositions = 0;
		long sharedPositions = 0;
		long comparedPositions = 0;
		long rawCopySegments = 0;
		long goodCopySegments = 0;

		try
		{
			// Declare output file //
	 	 	PrintStream outCopySegments = null; // declare a print stream object for copynumber segments

	 		outCopySegments = new PrintStream( new FileOutputStream(outputName + ".copynumber") );
	 		outCopySegments.println("chrom\tchr_start\tchr_stop\tnum_positions\tnormal_depth\ttumor_depth\tlog2_ratio\tgc_content");

	 		// Prepare file readers for normal and tumor pileups //

	 		BufferedReader normal = new BufferedReader(new FileReader(normalPileupFile));
		    BufferedReader tumor = new BufferedReader(new FileReader(tumorPileupFile));

		    if(!(normal.ready() && tumor.ready()))
		    {
		    	// Delay a few seconds to let SAMtools pileup start outputting //
		    	try {
			    	Thread.sleep(5000);

			    	if(!(normal.ready() && tumor.ready()))
			    		Thread.sleep(5000);

			    	if(!(normal.ready() && tumor.ready()))
			    		Thread.sleep(5000);

			    	if(!(normal.ready() && tumor.ready()))
			    		Thread.sleep(5000);
		    	}
		    	catch(Exception e)
		    	{

		    	}

		    }

		    // Exit if files not ready after waiting //

		    if(!(normal.ready() && tumor.ready()))
		    {
		    	System.err.println("ERROR: Invalid input file(s)");
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

		    // Parameters for copy number calling //
		    String copyChrom = "";
		    int copyStart = 0;
		    int copyStop = 0;
		    int copyDepthNormal = 0;
		    int copyDepthTumor = 0;
		    long copySumNormal = 0;
		    long copySumTumor = 0;
		    long copyPositions = 0;
		    long copyPositionsGC = 0;

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

					    	// If still less than tumor position, look for homozygous del //
					    	if(posNormal < posTumor)
					    	{
					    		int pileupDepthNormal = 0;
					    		String normalQualities = "";

			    				// Pileup Files have 6-7 columns //
		    					if(normalContents.length <= 7)
		    					{
		    						pileupDepthNormal = Integer.parseInt(normalContents[3]);
		    						normalQualities = normalContents[5];
		    					}
		    					// Pileup lines in CNS files have 10-11 columns
		    					else if (normalContents.length >= 10 && normalContents.length <= 11)
		    					{
		    						pileupDepthNormal = Integer.parseInt(normalContents[7]);
		    						normalQualities = normalContents[9];
		    					}

					    	}
					    	else
					    	{

					    	}
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
		    			refBase = tumorContents[2];

//		    			 Parse out base qualities //
	    				String[] normalContents = lineNormal.split("\t");
	    				int pileupDepthNormal = 0;
	    				int pileupDepthTumor = 0;
	    				String normalQualities = "";
	    				String tumorQualities = "";

	    				// Pileup Files have 6-7 columns //
    					if(normalContents.length >= 6 && normalContents.length <= 7)
    					{
    						pileupDepthNormal = Integer.parseInt(normalContents[3]);
    						normalQualities = normalContents[5];
    					}
    					// Pileup lines in CNS files have 10-11 columns
    					else if (normalContents.length >= 10 && normalContents.length <= 11)
    					{
    						pileupDepthNormal = Integer.parseInt(normalContents[7]);
    						normalQualities = normalContents[9];
    					}

	    				// Pileup Files have 6-7 columns //
    					if(tumorContents.length >= 6 && tumorContents.length <= 7)
    					{
    						tumorQualities = tumorContents[5];
    						pileupDepthTumor = Integer.parseInt(tumorContents[3]);
    					}
    					// Pileup lines in CNS files have 10-11 columns
    					else if (tumorContents.length >= 10 && tumorContents.length <= 11)
    					{
    						tumorQualities = tumorContents[9];
    						pileupDepthTumor = Integer.parseInt(tumorContents[7]);
    					}

    					// If either sample met the minimum coverage and both had at least one read //

//    					if((pileupDepthNormal >= minCoverage || pileupDepthTumor >= minCoverage) && normalQualities.length() > 0 && tumorQualities.length() > 0)

    					// We want the normal sample to meet the minimum coverage because that's the comparator //
    					if(pileupDepthNormal >= minCoverage && normalQualities.length() > 0) // && tumorQualities.length() > 0)
    					{
    						comparedPositions++;
//    						 Get the depth of bases above minimum quality //

    	    				int normalDepth = VarScan.qualityDepth(normalQualities, minBaseQual);
    	    				int tumorDepth = VarScan.qualityDepth(tumorQualities, minBaseQual);

    	    				// Determine if we have a copy changepoint //
    	    				// If this base is not contiguous with the copyRegion
    	    				// If the normal or tumor depth changes //

    	    				int diffNormal = Math.abs(copyDepthNormal - normalDepth);
    	    				int diffTumor = Math.abs(copyDepthTumor - tumorDepth);
    	    				int posDiff = posTumor - copyStop;

    	    				// DETERMINE IF WE CONTINUE THIS REGION OR PROCESS IT AND START A NEW ONE //

    	    				boolean continueFlag = false;

    	    				// If chromosomes differ or contiguity broken, process the region //

    	    				if(posDiff > 2 || !(copyChrom.equals(chromTumor)))
    	    				{
    	    					continueFlag = false;
    	    				}
    	    				else
    	    				{
    	    					if(copyPositions >= maxSegmentSize)
    	    					{
    	    						continueFlag = false;
    	    					}
    	    					else if(diffNormal <= 2 && diffTumor <= 2)
    	    					{
    	    						continueFlag = true;
    	    					}
    	    					else
    	    					{
    	    						// Do a Fisher's exact test on the copy number changes. ##

        	    					double changePvalue = VarScan.getSignificance(copyDepthNormal, copyDepthTumor, normalDepth, tumorDepth);

        	    					// If depth change not significant, continue with region //
        	    					if(changePvalue >= pValueThreshold)
        	    					{
        	    						continueFlag = true;
        	    					}
        	    					else
        	    					{
        	    						continueFlag = false;
        	    					}

    	    					}
    	    				}


    	    				// If continuing, extend this region and don't process yet //

    	    				if(continueFlag)
    	    				{
    	    					copySumNormal += normalDepth;
    	    					copySumTumor += tumorDepth;
    	    					copyPositions++;
    	    					if(refBase.equals("G") || refBase.equals("C") || refBase.equals("g") || refBase.equals("c"))
    	    						copyPositionsGC++;
    	    					copyStop = posTumor;
    	    				}

    	    				// Otherwise, process this region (if it qualifies) and start a new one //

    	    				else
    	    				{
    	    					if(copyPositions >= minSegmentSize)
    	    					{
    	    						rawCopySegments++;
    	    						String regionResults = processCopyRegion(copyChrom, copyStart, copyStop, copyPositions, copyPositionsGC, copySumNormal, copySumTumor, minCoverage, dataRatio);

    	    						if(regionResults.length() > 0)
    	    						{
    	    							outCopySegments.println(regionResults);
    	    							goodCopySegments++;
    	    						}
    	    					}

    	    					// Start a new copyNumber region //
    	    					copyChrom = chromTumor;
    	    					copyStart = posTumor;
    	    					copyStop = posTumor;
    	    					copyDepthNormal = normalDepth;
    	    					copyDepthTumor = tumorDepth;
    	    					copySumNormal = normalDepth;
    	    					copySumTumor = tumorDepth;
    	    					copyPositions = 1;
    	    					if(refBase.equals("G") || refBase.equals("C") || refBase.equals("g") || refBase.equals("c"))
    	    						copyPositionsGC = 1;
    	    					else
    	    						copyPositionsGC = 0;
    	    				}


    					}
    					else
    					{
    						// If minimum coverage was not met, print region //
	    					// If we had a copyNumber region that met minimum coverage, report it //
	    					if(copyPositions >= minSegmentSize)
	    					{
	    						rawCopySegments++;
	    						String regionResults = processCopyRegion(copyChrom, copyStart, copyStop, copyPositions, copyPositionsGC, copySumNormal, copySumTumor, minCoverage, dataRatio);

	    						if(regionResults.length() > 0)
	    						{
	    							outCopySegments.println(regionResults);
	    							goodCopySegments++;
	    						}
	    					}

	    					// Reset the copyNumber region //
	    					copyChrom = "";
	    					copyStart = 0;
	    					copyStop = 0;
	    					copyDepthNormal = 0;
	    					copyDepthTumor = 0;
	    					copySumNormal = 0;
	    					copySumTumor = 0;
	    					copyPositions = 0;
	    					copyPositionsGC = 0;
    					}

	    				// Record this chromosome //

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

			// If we had a copyNumber region that met minimum coverage, report it //
			if(copyPositions > minSegmentSize)
			{
				rawCopySegments++;
				String regionResults = processCopyRegion(copyChrom, copyStart, copyStop, copyPositions, copyPositionsGC, copySumNormal, copySumTumor, minCoverage, dataRatio);

				if(regionResults.length() > 0)
				{
					outCopySegments.println(regionResults);
					goodCopySegments++;
				}
			}


		    outCopySegments.close();

		    System.err.println(tumorPositions + " positions in tumor");
		    System.err.println(sharedPositions + " positions shared in normal"); //stats.get("sharedPositions")
		    System.err.println(comparedPositions + " had sufficient coverage for comparison"); //stats.get("comparedPositions")

		    System.err.println(rawCopySegments + " raw copynumber segments with size > " + minSegmentSize);
		    System.err.println(goodCopySegments + " good copynumber segments with depth > " + minCoverage);
		}
		catch (IOException e)
		{
			System.err.println("File Parsing Exception: " + e.getLocalizedMessage());
			e.printStackTrace(System.err);
			System.exit(11);
		}
	}


	/**
	 * Calculates relative tumor copynumber for a contiguous segment
	 *
	 * @param	args	Command-line arguments
	 * @return			HashMap of parameter names and their values
	 */
	static String processCopyRegion(String copyChrom, int copyStart, int copyStop, long copyPositions, long copyPositionsGC, long copySumNormal, long copySumTumor, int minCoverage, double dataRatio)
	{
	    DecimalFormat oneDigit = new DecimalFormat("#0.0");
	    DecimalFormat threeDigits = new DecimalFormat("#0.000");

		try
		{
//			 Calculate average depth //
			float avgNormal = (float) copySumNormal / (float) copyPositions;
			float avgTumor = (float) copySumTumor / (float) copyPositions;
			// Adjust tumor depth for ratio
			float adjustedTumorDepth = (float) dataRatio * (float) avgTumor;

			float gcContent = (float) copyPositionsGC / (float) copyPositions * 100;

			if(avgNormal >= minCoverage || avgTumor >= minCoverage)
			{
    			// Determine ratio and diff //
    			if(avgNormal >= 0.01 && avgTumor >= 0.01)
    			{
    		 		float tumorNormalRatio = adjustedTumorDepth / avgNormal;
    				double log2ratio = Math.log(tumorNormalRatio) / Math.log(2);

					return(copyChrom + "\t" + copyStart + "\t" + copyStop + "\t" + copyPositions + "\t" + oneDigit.format(avgNormal) + "\t" + oneDigit.format(avgTumor) + "\t" + threeDigits.format(log2ratio) + "\t" + oneDigit.format(gcContent));
    			}
    			else if (avgTumor >= 0.01)
    			{
    				// If only tumor has coverage, handle it //
    				double log2ratio = 2.00;
					return(copyChrom + "\t" + copyStart + "\t" + copyStop + "\t" + copyPositions + "\t" + oneDigit.format(avgNormal) + "\t" + oneDigit.format(avgTumor) + "\t" + threeDigits.format(log2ratio) + "\t" + oneDigit.format(gcContent));
    			}
    			else
    			{
    				// If only normal has coverage, mark as homozygyous deletion //
    				double log2ratio = -2.00;
					return(copyChrom + "\t" + copyStart + "\t" + copyStop + "\t" + copyPositions + "\t" + oneDigit.format(avgNormal) + "\t" + oneDigit.format(avgTumor) + "\t" + threeDigits.format(log2ratio) + "\t" + oneDigit.format(gcContent));
    			}

			}
			else
			{
//				System.err.println("Warning: Not reporting region " + copyChrom + " " + copyStart + " " + copyStop + " " + copyPositions + " " + avgNormal + " " + avgTumor);
			}
		}
		catch(Exception e)
		{
			System.err.println("Warning: Error while processing copynumber segment:" + e.getMessage());
		}


		return("");
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
