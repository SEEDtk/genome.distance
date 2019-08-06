/**
 *
 */
package org.theseed.genome.distance;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;
import org.theseed.genome.Genome;
import org.theseed.genome.GenomeDirectory;
import org.theseed.proteins.RoleMap;

/**
 * Compare a base genome to all the genomes in a directory.  The output contains
 * the percent similarity for each input genome by genome ID.  The positional parameters are
 * the name of the base GTO file and the name of the genome directory.  Multiple
 * genome directories can be specified.
 *
 * @author Bruce Parrello
 *
 */
public class DistanceProcessor {

	// FIELDS
	/** map of useful roles */
	RoleMap usefulRoles;
	/** current input genome directory */
	GenomeDirectory gtoDir;
	/** base genome */
	Genome baseGenome;

    // COMMAND LINE

    /** help option */
    @Option(name="-h", aliases={"--help"}, help=true)
    private boolean help;

    /** TRUE if we want progress messages */
    @Option(name="-v", aliases={"--verbose", "--debug"}, usage="display progress on STDERR")
    private boolean debug;

    /** name of the useful-role file */
    @Option(name="-R", aliases={"--roles"}, metaVar="roleFile", required=true,
    		usage="file of useful roles")
    private File roleFile;

    /** name of the base genome's GTO file */
    @Argument(index=0, metaVar="gtoFile", required=true, usage="base genome GTO file")
    private File gtoFile;

    /** name of the input genome directory */
    @Argument(index=1, metaVar="gtoDir1 gtoDir2 ...", required=true, usage="directory of input GTOs")
    private List<File> genomeDirs;

    /**
     * Parse command-line options to specify the parameters of this object.
     *
     * @param args	an array of the command-line parameters and options
     *
     * @return TRUE if successful, FALSE if the parameters are invalid
     */
    public boolean parseCommand(String[] args) {
        boolean retVal = false;
        // Set the defaults.
        this.help = false;
        this.debug = false;
        // Parse the command line.
        CmdLineParser parser = new CmdLineParser(this);
        try {
            parser.parseArgument(args);
            if (this.help) {
                parser.printUsage(System.err);
            } else {
            	// Get the role file.
            	this.usefulRoles = RoleMap.load(this.roleFile);
            	// Get the base genome.
            	this.baseGenome = new Genome(this.gtoFile);
            	// Verify the genome directories.
            	for (File genomeDir : this.genomeDirs) {
            		if (! genomeDir.isDirectory()) {
            			throw new IOException("Directory " + genomeDir + " not found or invalid.");
            		}
            	}
                // We made it this far, we can run the application.
                retVal = true;
            }
        } catch (CmdLineException e) {
            System.err.println(e.getMessage());
            // For parameter errors, we display the command usage.
            parser.printUsage(System.err);
        } catch (IOException e) {
            System.err.println(e.getMessage());
        }
        return retVal;
    }

    public void run() {
    	try {
	    	// Write the output header.
	    	System.out.println("genome_id\tgenome_name\tsimilarity");
	    	// Create the main measurement object.
	    	if (debug) System.err.println("Loading genome " + this.baseGenome + ".");
	    	Measurer baseKmers = new Measurer(this.baseGenome, this.usefulRoles);
	    	// Loop through the input directories.
	    	for (File inDir : this.genomeDirs) {
	    		if (this.debug) System.err.println("Processing directory " + inDir + ".");
	    		GenomeDirectory genomes = new GenomeDirectory(inDir);
	    		// Loop through the genomes.
	    		for (Genome genome : genomes) {
	    			if (this.debug) System.err.println("Processing genome " + genome + ".");
	    			double percent = baseKmers.computePercentSimilarity(genome);
	    			System.out.format("%s\t%s\t%8.2f%n", genome.getId(),
	    					genome.getName(), percent);
	    		}
	    	}
    	} catch (IOException e) {
    		System.err.println(e.getMessage());
    	}
    }
}
