package net.sf.varscan;

import java.io.BufferedReader;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.PrintStream;
import java.util.BitSet;
import java.util.HashMap;

/**
 * A class for processing VarScan output by somatic status and confidence
 *
 * @version	2.3
 *
 * @author Daniel C. Koboldt <dkoboldt@genome.wustl.edu>
 *
 */
public class ProcessSomatic {

	public ProcessSomatic(String[] args)
	{
		String usage = "USAGE: java -jar VarScan.jar process [status-file] OPTIONS\n" +
		"\tstatus-file - The VarScan output file for SNPs or Indels\n" +
		"\tOPTIONS\n" +
		"\t--min-tumor-freq - Minimum variant allele frequency in tumor [0.10]\n" +
		"\t--max-normal-freq - Maximum variant allele frequency in normal [0.05]\n" +
		"\t--p-value - P-value for high-confidence calling [0.07]";

		// Parse command-line parameters //
		HashMap<String, String> params = VarScan.getParams(args);

		// Try adjusting any provided parameters based on user inut //
		String statusFile = "varScan.output";

		double maxNormalFreq = 0.05;
		double minTumorFreq = 0.10;
		double pValueForHC = 0.07;

		try
		{
			if(params.containsKey("min-tumor-freq"))
				 minTumorFreq = Double.parseDouble(params.get("min-tumor-freq"));
			if(params.containsKey("max-normal-freq"))
				 maxNormalFreq = Double.parseDouble(params.get("max-normal-freq"));
			if(params.containsKey("p-value"))
				 pValueForHC = Double.parseDouble(params.get("p-value"));
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

		try
		{
			//			 Obtain input file //
			BufferedReader infile = VarScan.getInfile(args);

	    	// If no input, print usage //

	    	if(infile == null)
	    	{
	    		System.out.println(usage);
				return;
	    	}

	    	// If input file not ready, give it a few seconds //
	    	int numNaps = 0;

	    	while(!infile.ready())
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

			if(args.length > 1)
				statusFile = args[1];


			// Declare output files //
			PrintStream outSomatic = null;
			PrintStream outSomaticHC = null;
			PrintStream outGermline = null;
			PrintStream outGermlineHC = null;
			PrintStream outLOH = null;
			PrintStream outLOHHC = null;

			boolean isVCF = false;

			// Open output files for Somatic, Somatic HC, Germline, and LOH //

			if(statusFile.endsWith(".vcf"))
			{
				isVCF = true;
				String nameString = statusFile.replace(".vcf", "");
				System.err.println("Opening output files: " + nameString + ".Somatic.vcf " + nameString + ".Germline.vcf " + nameString + ".LOH.vcf ");
				outSomatic = new PrintStream( new FileOutputStream(nameString + ".Somatic.vcf") );
				outSomaticHC = new PrintStream( new FileOutputStream(nameString + ".Somatic.hc.vcf") );
				outGermline = new PrintStream( new FileOutputStream(nameString + ".Germline.vcf") );
				outGermlineHC = new PrintStream( new FileOutputStream(nameString + ".Germline.hc.vcf") );
				outLOH = new PrintStream( new FileOutputStream(nameString + ".LOH.vcf") );
				outLOHHC = new PrintStream( new FileOutputStream(nameString + ".LOH.hc.vcf") );

			}
			else
			{
				System.err.println("Opening output files: " + statusFile + ".Somatic " + statusFile + ".Germline " + statusFile + ".LOH ");
				outSomatic = new PrintStream( new FileOutputStream(statusFile + ".Somatic") );
				outSomaticHC = new PrintStream( new FileOutputStream(statusFile + ".Somatic.hc") );
				outGermline = new PrintStream( new FileOutputStream(statusFile + ".Germline") );
				outGermlineHC = new PrintStream( new FileOutputStream(statusFile + ".Germline.hc") );
				outLOH = new PrintStream( new FileOutputStream(statusFile + ".LOH") );
				outLOHHC = new PrintStream( new FileOutputStream(statusFile + ".LOH.hc") );

			}





			// Reset counters //
			int numProcessed = 0;
			int numSomatic = 0;
			int numSomaticHC = 0;
			int numGermline = 0;
			int numLOH = 0;
			int numGermlineHC = 0;
			int numLOHHC = 0;

	 	 	// Declare file-parsing variables //
	 	 	String line = "";
	 	 	int lineCounter = 0;

	 	 	// Actually parse the infile //

	 	 	while ((line = infile.readLine()) != null)
		    {
		    	lineCounter++;

	    		String[] lineContents = line.split("\t");

	    		if(lineContents.length >= 1)
	    		{
	    			// Try to parse chrom and position //
	    			try
	    			{
	    				String refName = lineContents[0];

	    				if(line.startsWith("#"))
		    			{
		    				// VCF header line //
		    				isVCF = true;
		    				outSomatic.println(line);
		    				outSomaticHC.println(line);
		    				outGermline.println(line);
		    				outGermlineHC.println(line);
		    				outLOH.println(line);
		    				outLOHHC.println(line);
		    			}
	    				else if(refName.equals("chrom") || refName.equals("Chrom"))
	    				{
		    				outSomatic.println(line);
		    				outSomaticHC.println(line);
		    				outGermline.println(line);
		    				outGermlineHC.println(line);
		    				outLOH.println(line);
		    				outLOHHC.println(line);
	    				}
	    				else
	    				{

		    				int position = Integer.parseInt(lineContents[1]);
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
	    					else if(refName.equals("chrom"))
	    					{
	    						// VarScan header //
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

		    				if(normalReads1 > 0 || normalReads2 > 0)
		    					normalFreq = (double) normalReads2 / (double) (normalReads1 + normalReads2);

		    				if(tumorReads1 > 0 || tumorReads2 > 0)
		    					tumorFreq = (double) tumorReads2 / (double) (tumorReads1 + tumorReads2);

		    				numProcessed++;

		    				if(somaticStatus.equals("Somatic") || somaticStatus.equals("2"))
		    				{
		    					numSomatic++;
		    					outSomatic.println(line);
		    					if(normalFreq <= maxNormalFreq && tumorFreq >= minTumorFreq && somaticPvalue <= pValueForHC)
		    					{
		    						numSomaticHC++;
		    						outSomaticHC.println(line);
		    					}
		    				}
		    				else if(somaticStatus.equals("Germline") || somaticStatus.equals("1"))
		    				{
		    					numGermline++;
		    					outGermline.println(line);
		    					if(normalFreq >= minTumorFreq && tumorFreq >= minTumorFreq && somaticPvalue <= pValueForHC)
		    					{
		    						numGermlineHC++;
		    						outGermlineHC.println(line);
		    					}
		    				}
		    				else if(somaticStatus.equals("LOH") || somaticStatus.equals("3"))
		    				{
		    					numLOH++;
		    					outLOH.println(line);
		    					double normalHetDistance = Math.abs(0.50 - normalFreq);
		    					double tumorHetDistance = Math.abs(0.50 - tumorFreq);
		    					if(normalFreq >= minTumorFreq && tumorHetDistance > normalHetDistance && somaticPvalue <= pValueForHC)
		    					{
		    						numLOHHC++;
		    						outLOHHC.println(line);
		    					}
		    				}
	    				}




	    			}
	    			catch(Exception e)
	    			{
	    				if(lineCounter == 1)
	    				{
	    					// Print the header to the output files //
	    					outSomatic.println(line);
	    					outSomaticHC.println(line);
	    					outGermline.println(line);
	    					outLOH.println(line);

	    				}
	    				else
	    					System.err.println("Warning: Unable to parse chrom/position from " + line);
	    			}
	    		}
		    }

	 	 	// Close the files //
	 	 	infile.close();
	 	 	outSomatic.close();
	 	 	outSomaticHC.close();
	 	 	outGermline.close();
	 	 	outLOH.close();
	 	 	outGermlineHC.close();
	 	 	outLOHHC.close();

			// Print the status report //

			System.out.println(numProcessed + " VarScan calls processed");
			System.out.println(numSomatic + " were Somatic (" + numSomaticHC + " high confidence)");
			System.out.println(numGermline + " were Germline (" + numGermlineHC + " high confidence)");
			System.out.println(numLOH + " were LOH (" + numLOHHC + " high confidence)");

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