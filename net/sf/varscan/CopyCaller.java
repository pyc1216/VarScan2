/**
 * @(#)CopyCaller.java
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
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.BitSet;
import java.util.HashMap;

/**
 * A class for calling/GC-adjusting copy number variants from raw somatic copynumber output
 *
 * @version	2.3
 *
 * @author Daniel C. Koboldt <dkoboldt@genome.wustl.edu>
 *
 */
public class CopyCaller {
	public CopyCaller(String[] args, HashMap<String, String> params)
	{
		String usage = "USAGE: java -jar VarScan.jar copyCaller [varScan.copynumber] OPTIONS\n" +
		"This command will adjust VarScan copynumber output for GC content, apply amp/del thresholds,\n and (optionally) recenter the data\n" +
		"\tINPUT:\n" +
		"\tRaw output from the VarScan copynumber command (eg. varScan.output.copynumber)\n\n" +
		"\tOPTIONS:\n" +
		"\t--output-file\tOutput file to contain the calls\n" +
		"\t--output-homdel-file\tOptional output file for candidate homozygous deletions\n" +
		"\t--min-coverage\tMinimum normal read depth at a position to make a call [20]\n" +
		"\t--min-tumor-coverage\tMinimum tumor read depth at a position to make a non-homdel call [10]\n" +
		"\t--max-homdel-coverage\tMaximum depth in tumor for candidate homozygous deletions [5]\n" +
		"\t--amp-threshold\tLower bound for log ratio to call amplification [0.25]\n" +
		"\t--del-threshold\tUpper bound for log ratio to call deletion (provide as positive number) [0.25]\n" +
		"\t--min-region-size\tMinimum size (in bases) for a region to be counted [10]\n" +
		"\t--recenter-up\tRecenter data around an adjusted baseline > 0 [0]\n" +
		"\t--recenter-down\tRecenter data around an adjusted baseline < 0 [0]\n";

		// Set parameter defaults //

		String regionsFile = "";
		String outputFile = "";
		String homdelFile = "";
		int minCoverage = 20;
		int minTumorCoverage = 10;
		int maxHomdelCoverage = 5;
		int minRegionSize = 10;
		double ampThreshold = 0.25;
		double delThreshold = -0.25;
		double recenterBaseline = 0.00;
		Float[] gcLogSum = new Float[101];
		Integer[] gcLogNum = new Integer[101];

		// Reset GC bin //

		for(int i = 0; i <= 100; i++)
		{
			gcLogSum[i] = (float) 0;
			gcLogNum[i] = 0;
		}

		// Get any user-provided parameters //

		try
		{
			if(params.containsKey("min-coverage"))
				 minCoverage = Integer.parseInt(params.get("min-coverage"));

			if(params.containsKey("min-tumor-coverage"))
				 minTumorCoverage = Integer.parseInt(params.get("min-tumor-coverage"));

			if(params.containsKey("max-homdel-coverage"))
				 maxHomdelCoverage = Integer.parseInt(params.get("max-homdel-coverage"));

			if(params.containsKey("min-region-size"))
				 minRegionSize = Integer.parseInt(params.get("min-region-size"));

			if(params.containsKey("amp-threshold"))
				ampThreshold = Double.parseDouble(params.get("amp-threshold"));

			if(params.containsKey("del-threshold"))
			{
				delThreshold = 0 - Double.parseDouble(params.get("del-threshold"));
			}

			if(params.containsKey("recenter-up"))
			{
				recenterBaseline = Double.parseDouble(params.get("recenter-up"));
			}

			if(params.containsKey("recenter-down"))
			{
				recenterBaseline = Double.parseDouble(params.get("recenter-down"));
				recenterBaseline = 0.00 - recenterBaseline;
			}

			if(params.containsKey("output-file"))
				 outputFile = params.get("output-file");

			if(params.containsKey("output-homdel-file"))
				 homdelFile = params.get("output-homdel-file");

			 System.err.println("Min coverage:\t" + minCoverage);
		}
		catch(Exception e)
		{
	    	System.err.println("Input Parameter Threw Exception: " + e.getLocalizedMessage());
	    	System.err.println("Parsing " + params.get("del-threshold"));
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
	    stats.put("numRegions", 0);
	    stats.put("metMinDepth", 0);
	    stats.put("metMinSize", 0);
	    stats.put("numAmp", 0);
	    stats.put("numDel", 0);
	    stats.put("numHomDel", 0);
	    stats.put("numNeutral", 0);

	    DecimalFormat threeDigits = new DecimalFormat("#0.000");

	    HashMap<String, Long> baseCounts = new HashMap<String, Long>();
	    baseCounts.put("numAmp", (long) 0);
	    baseCounts.put("numDel", (long) 0);
	    baseCounts.put("numHomDel", (long) 0);
	    baseCounts.put("numNeutral", (long) 0);

	    try
	    {
	    	// If output file was provided, open it //
	    	PrintStream out = null; // declare a print stream object
	    	PrintStream outGC = null;
	    	PrintStream outHomdel = null; // declare a print stream object

	    	if(params.containsKey("output-file"))
	    	{
	    		out = new PrintStream( new FileOutputStream(outputFile) );
    			out.println("chrom\tchr_start\tchr_stop\tnum_positions\tnormal_depth\ttumor_depth\tadjusted_log_ratio\tgc_content\tregion_call\traw_ratio");

    			outGC = new PrintStream ( new FileOutputStream(outputFile + ".gc") );
    			outGC.println("gc\tregions\tavg_log2\tmean_sd_log2");
	    	}

	    	if(params.containsKey("output-homdel-file"))
	    	{
	    		outHomdel = new PrintStream( new FileOutputStream(homdelFile) );
    			outHomdel.println("chrom\tchr_start\tchr_stop\tnum_positions\tnormal_depth\ttumor_depth\tadjusted_log_ratio\tgc_content\tregion_call\traw_ratio");
	    	}

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


			boolean gcWarned = false;
	    	// If input input was provided, begin parsing it //

	    	if(in.ready())
	    	{
	    		while ((line = in.readLine()) != null)
	    		{
	    			// Output progress line //
	    			if(params.containsKey("verbose") && (stats.get("numRegions") % 10000) == 0)
		        		System.err.println(stats.get("numRegions") + " regions parsed...");

	    			// Begin try-catch for line parsing //

	    			try
	    			{
	    				String[] lineContents = line.split("\t");

	    				// Verify expected copynumber regions format //

	    				if(lineContents.length > 4 && lineContents[0].length() > 0 && lineContents[1].length() > 0 && lineContents[2].length() > 0 && lineContents[3].length() > 0)
	    				{
	    					String refName = lineContents[0];

	    					// Print the header line //
	    					if(refName.equals("chrom"))
	    					{

	    					}
	    					else if(!refName.equals("chrom"))
	    					{
	    						stats.put("numRegions", (stats.get("numRegions") + 1));

		    					long regionStart = Long.parseLong(lineContents[1]);
		    					long regionStop = Long.parseLong(lineContents[2]);
		    					long numPositions = Long.parseLong(lineContents[3]);

		    					// Fix locale-parsing issues //
		    					float normalDepth = Float.parseFloat(lineContents[4].replace(',', '.'));
		    					float tumorDepth = Float.parseFloat(lineContents[5].replace(',', '.'));
			    	        	double logratio = Double.parseDouble(lineContents[6].replace(',', '.'));

			    	        	if(recenterBaseline != 0)
			    	        		logratio = logratio - recenterBaseline;

			    	        	if(lineContents.length >= 8)
			    	        	{
			    	        		// Apply coverage threshold //
			    	        		if(normalDepth >= minCoverage && tumorDepth >= minTumorCoverage)
			    	        		{
				    	        		float gcContent = Float.parseFloat(lineContents[7].replace(',', '.'));
				    	        		int gcBin = (int) gcContent;
				    	        		if(gcBin >= 0 && gcBin <= 100)
				    	        		{
				    	        			gcLogSum[gcBin] += (float) logratio;
				    	        			gcLogNum[gcBin]++;
				    	        		}
			    	        		}

			    	        	}
			    	        	else
			    	        	{
			    	        		// No GC information.. warn if we haven't already //
			    	        		if(!gcWarned)
			    	        		{
			    	        			System.err.println("Warning: Older VarScan copynumber output (without GC content column) detected, so no GC adjustment will be performed");
			    	        			gcWarned = true;
			    	        		}
			    	        	}


	    					}
	    					else
	    					{
	    						// Don't process the header line //
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
		    				System.exit(11);
	    				}

	    		    }


	    		}

	    	}
	    	else
	    	{
	    		System.err.println("Input was not ready for parsing!");
		    	System.exit(10);
	    	}

	    	in.close();

	    	// Get overall mean copy number //
	    	float totalAvgSum = (float) 0;
	    	long totalAvgNum = 0;
	    	// Print GC content by bin //
			for(int i = 0; i <= 100; i++)
			{
				if(gcLogNum[i] > 0)
				{
					totalAvgSum += gcLogSum[i];
					totalAvgNum += gcLogNum[i];
				}
			}

			float totalAvgLog = totalAvgSum / (float) totalAvgNum;

			if(!gcWarned)
			{
				System.err.println(totalAvgLog + " was the average log2 of copy number change");
				System.err.println("Copy number change by GC content bin:");
				System.err.println("bin\tregions\tavg_log2\tmean_sd");
			}

			Float[] gcLogMeanSD = new Float[101];
			try
			{
		    	// Print GC content by bin //
				for(int i = 0; i <= 100; i++)
				{
					if(gcLogNum[i] > 0)
					{
						float binAvgLog = gcLogSum[i] / (float) gcLogNum[i];
						float binMeanSD = binAvgLog - totalAvgLog;
						outGC.println(i + "\t" + gcLogNum[i] + "\t" + binAvgLog + "\t" + binMeanSD);
						System.err.println(i + "\t" + gcLogNum[i] + "\t" + binAvgLog + "\t" + binMeanSD);
						gcLogMeanSD[i] = binMeanSD;
					}
					else
					{
						gcLogMeanSD[i] = (float) 0;
					}
				}
			}
			catch(Exception e)
			{
				System.err.println("Insufficient data for GC adjustment, so adjustment will be skipped");
				for(int i = 0; i <= 100; i++)
				{
					gcLogMeanSD[i] = (float) 0;
				}
			}



			// Re-parse the input file //

			in = VarScan.getInfile(args);
	    	// If input input was provided, begin parsing it //


	    	if(in.ready())
	    	{
	    		while ((line = in.readLine()) != null)
	    		{
	    			// Output progress line //
	    			if(params.containsKey("verbose") && (stats.get("numRegions") % 10000) == 0)
		        		System.err.println(stats.get("numRegions") + " regions parsed...");

	    			// Begin try-catch for line parsing //

	    			try
	    			{
	    				String[] lineContents = line.split("\t");

	    				// Verify expected copynumber regions format //

	    				if(lineContents.length > 4 && lineContents[0].length() > 0 && lineContents[1].length() > 0 && lineContents[2].length() > 0 && lineContents[3].length() > 0)
	    				{
	    					String refName = lineContents[0];

	    					// Print the header line //
	    					if(refName.equals("chrom"))
	    					{

	    					}
	    					else if(!refName.equals("chrom"))
	    					{
	    						stats.put("numRegions", (stats.get("numRegions") + 1));

		    					long regionStart = Long.parseLong(lineContents[1]);
		    					long regionStop = Long.parseLong(lineContents[2]);
		    					long numPositions = Long.parseLong(lineContents[3]);
		    					float normalDepth = Float.parseFloat(lineContents[4].replace(',', '.'));
		    					float tumorDepth = Float.parseFloat(lineContents[5].replace(',', '.'));
			    	        	double logratio = Double.parseDouble(lineContents[6].replace(',', '.'));
			    	        	double adjustedRatio = logratio;

			    	        	// If recentering, adjust the adjusted log ratio //

			    	        	if(recenterBaseline != 0)
			    	        		adjustedRatio = adjustedRatio - recenterBaseline;

			    	        	float gcContent = (float) -1;
			    	        	if(lineContents.length >= 8)
			    	        	{
			    	        		gcContent = Float.parseFloat(lineContents[7].replace(',', '.'));
			    	        		int gcBin = (int) gcContent;
			    	        		// If there was an adjustment for this GC bin, make it so //
			    	        		if(gcBin >= 0 && gcBin <= 100) // && normalDepth >= minCoverage && tumorDepth >= minTumorCoverage
			    	        		{
			    	        			if(gcLogMeanSD[gcBin] != (float) 0)
			    	        			{
			    	        				adjustedRatio = adjustedRatio - gcLogMeanSD[gcBin];
			    	        			}
			    	        		}
			    	        	}


			    	        	// Check to see if this position meets minimum depth //
		    	        		long regionSize = regionStop - regionStart + 1;

			    	        	if(normalDepth >= minCoverage && tumorDepth >= minTumorCoverage)
			    	        	{
			    	        		stats.put("metMinDepth", (stats.get("metMinDepth") + 1));

			    	        		String regionCall = "neutral";

			    	        		if(regionSize >= minRegionSize)
			    	        		{
			    	        			stats.put("metMinSize", (stats.get("metMinSize") + 1));

			    	        			// Determine class based on user-specified thresholds //

			    	        			if(adjustedRatio >= ampThreshold)
			    	        			{
			    	        				stats.put("numAmp", (stats.get("numAmp") + 1));
			    	        				baseCounts.put("numAmp", (baseCounts.get("numAmp") + regionSize));

			    	        				regionCall = "amp";
			    	        			}
			    	        			else if(adjustedRatio <= delThreshold)
			    	        			{
			    	        				stats.put("numDel", (stats.get("numDel") + 1));
			    	        				baseCounts.put("numDel", (baseCounts.get("numDel") + regionSize));
			    	        				regionCall = "del";
			    	        			}
			    	        			else
			    	        			{
			    	        				stats.put("numNeutral", (stats.get("numNeutral") + 1));
			    	        				baseCounts.put("numNeutral", (baseCounts.get("numNeutral") + regionSize));
			    	        			}

			    	        			String outLine = refName + "\t" + regionStart + "\t" + regionStop + "\t" + numPositions + "\t";
			    	        			outLine += normalDepth + "\t" + tumorDepth + "\t" + threeDigits.format(adjustedRatio) + "\t" + gcContent + "\t" + regionCall + "\t" + logratio;

			    	        			// Print to outfile or standardout //

			    	        			if(params.containsKey("output-file"))
			    	        			{
			    	        				out.println(outLine);
			    	        			}
			    	        			else
			    	        			{
			    	        				System.err.println(outLine);
			    	        			}
			    	        		}

			    	        	}
			    	        	else if(normalDepth >= minCoverage && tumorDepth <= maxHomdelCoverage && regionSize >= minRegionSize && adjustedRatio <= delThreshold)
			    	        	{
			    	        		// Output candidate homozygous deletion //
		    	        			String outLine = refName + "\t" + regionStart + "\t" + regionStop + "\t" + numPositions + "\t";
		    	        			outLine += normalDepth + "\t" + tumorDepth + "\t" + threeDigits.format(adjustedRatio) + "\t" + gcContent + "\thomozygous_deletion\t" + logratio;
		    	        			stats.put("numHomDel", (stats.get("numHomDel") + 1));
	    	        				baseCounts.put("numHomDel", (baseCounts.get("numHomDel") + regionSize));
	    	        				if(params.containsKey("output-homdel-file"))
	    	        				{
	    	        					outHomdel.println(outLine);
	    	        				}

			    	        	}

	    					}
	    					else
	    					{
	    						// Don't process the header line //
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
	    		    	System.err.println("Parsing Exception on line:\n" + line + "\n" + e.getMessage() + "\n" + e.getLocalizedMessage());
	    		    	System.err.println(e.toString());
	    				stats.put("numParsingExceptions", (stats.get("numParsingExceptions") + 1));

	    				if(stats.get("numParsingExceptions") >= 5)
	    				{
	    					System.err.println("Too many parsing exceptions encountered; exiting");
		    				System.exit(11);
	    				}

	    		    }


	    		}

	    	}
	    	else
	    	{
	    		System.err.println("Input was not ready for parsing!");
		    	System.exit(10);
	    	}

	    	if(params.containsKey("output-file"))
	    	{
	    		out.close();
	    		outGC.close();
	    	}



	    	// Print summary statistics //
	    	System.err.println(stats.get("numRegions") + " raw regions parsed");
	    	System.err.println(stats.get("metMinDepth") + " met min depth");
	    	System.err.println(stats.get("metMinSize") + " met min size");
	    	System.err.println(stats.get("numAmp") + " regions (" + baseCounts.get("numAmp") + " bp)" + " were called amplification (log2 > " + ampThreshold + ")");
	    	System.err.println(stats.get("numNeutral") + " regions (" + baseCounts.get("numNeutral") + " bp)" + " were called neutral");
	    	System.err.println(stats.get("numDel") + " regions (" + baseCounts.get("numDel") + " bp)" + " were called deletion (log2 <" + delThreshold + ")");
	    	System.err.println(stats.get("numHomDel") + " regions (" + baseCounts.get("numHomDel") + " bp)" + " were called homozygous deletion (normal cov >= " + minCoverage + " and tumor cov <= " + maxHomdelCoverage + ")");
	    }
	    catch(Exception e)
	    {
	    	System.err.println("Error parsing input: " + e.getLocalizedMessage());
	    	e.printStackTrace(System.err);
	    	System.exit(11);
	    }
	}


/**
 * Processes a region of copy number calls
 *
 * @param	regionRef	Reference or chromosome name
 * @param	regionStart	Start position on reference
 * @param	regionStop	Stop position on reference
 */
	static String processRegion(String regionRef, long regionStart, long regionStop, int regionCalls, long regionSumNormal, long regionSumTumor)
	{
		// Calculate average region depth //

		float regionDepthNormal = regionSumNormal / regionCalls;
		float regionDepthTumor = regionSumTumor / regionCalls;

		double tumorNormalRatio = 100;
		double log2ratio = 0.00;
		if(regionDepthNormal > 0)
		{
			tumorNormalRatio = regionDepthTumor / regionDepthNormal;
			log2ratio = Math.log(tumorNormalRatio) / Math.log(2);
		}

		return(regionRef + "\t" + regionStart + "\t" + regionStop + "\t" + regionCalls + "\t" + regionDepthNormal + "\t" + regionDepthTumor + "\t" + log2ratio);

	}
}