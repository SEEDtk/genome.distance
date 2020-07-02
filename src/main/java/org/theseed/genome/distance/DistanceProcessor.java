/**
 *
 */
package org.theseed.genome.distance;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.List;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Genome;
import org.theseed.genome.GenomeDirectory;
import org.theseed.proteins.RoleMap;
import org.theseed.proteins.RoleScanner;
import org.theseed.utils.BaseProcessor;

/**
 * Compare base genomes to all the genomes in a directory.  The output contains
 * the percent similarity for each input genome by genome ID.  The positional parameters are
 * the name of the directory containing the base GTO files and the name of the genome directory
 * containing the comparison genomes.  Multiple comparison genome directories can be specified.
 *
 * The command-line options are as follows.
 *
 * -v	show more detail on the log
 * -h	display usage information
 * -R	the name of a file containing the useful roles for the comparison, in the form of a role map
 * 		(tab delimited, role ID, comment, role description); if this parameter is omitted, the input
 * 		directory is scanned to create the map of useful roles
 *
 * --seed	ID of the seed protein; if this is specified, -R is required.  If specified, the seed protein
 * 			distances will be shown in a separate column
 *
 * @author Bruce Parrello
 *
 */
public class DistanceProcessor extends BaseProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(DistanceProcessor.class);

    /** map of useful roles */
    private RoleMap usefulRoles;

    // COMMAND LINE

    /** name of the useful-role file */
    @Option(name="-R", aliases={"--roles"}, metaVar="roleFile",
            usage="file of useful roles")
    private File roleFile;

    /** ID of the seed protein */
    @Option(name="--seed", metaVar="PhenTrnaSyntAlph", usage="ID of the seed protein")
    private String seedProtein;



    /** name of the base genome directory */
    @Argument(index=0, metaVar="gtoDir", required=true, usage="base genome GTO directory")
    private File baseDir;

    /** name of the input genome directory */
    @Argument(index=1, metaVar="gtoDir1 gtoDir2 ...", required=true, usage="directory of input GTOs")
    private List<File> genomeDirs;

    @Override
    protected void setDefaults() {
        this.seedProtein = null;
        this.roleFile = null;
    }

    @Override
    protected boolean validateParms() throws IOException {
        // Verify the base genome directory.
        if (! this.baseDir.isDirectory())
            throw new IOException("Directory " + this.baseDir + " not found or invalid.");
        // Verify the genome directories.
        for (File genomeDir : this.genomeDirs) {
            if (! genomeDir.isDirectory()) {
                throw new IOException("Directory " + genomeDir + " not found or invalid.");
            }
        }
        // Verify the role file.
        if (this.roleFile != null && ! this.roleFile.canRead())
            throw new FileNotFoundException("Role file " + this.roleFile + " not found or unreadable.");
        // Insure there is a role file if there is a seed protein.
        if (this.seedProtein != null && this.roleFile == null)
            throw new IOException("Seed protein requires specification of a role file.");
        // We made it this far, we can run the application.
        return true;
    }

    public void runCommand() throws Exception {
        // Get the base directory.
        GenomeDirectory baseGenomes = new GenomeDirectory(this.baseDir);
        // Create the role map.
        if (this.roleFile != null) {
            // Here we are reading in a role file.
            log.info("Reading roles from {}.", this.roleFile);
            this.usefulRoles = RoleMap.load(this.roleFile);
        } else {
            // Here we have to scan the input directory.
            log.info("Scanning genomes in {} to compute useful roles.", this.baseDir);
            RoleScanner roleMap = new RoleScanner();
            roleMap.addGenomes(baseGenomes);
            this.usefulRoles = roleMap;
        }
        // Write the output header.
        String header = "base_id\tbase_name\tgenome_id\tgenome_name\tsimilarity";
        if (this.seedProtein != null)
            header += "\tseed_distance";
        System.out.println(header);
        // Create the main measurement object.
        for (Genome baseGenome : baseGenomes) {
            log.info("Loading genome {}.", baseGenome);
            Measurer baseKmers = new Measurer(baseGenome, this.usefulRoles);
            // Loop through the input directories.
            for (File inDir : this.genomeDirs) {
                log.info("Processing directory {}.", inDir);
                GenomeDirectory genomes = new GenomeDirectory(inDir);
                // Loop through the genomes.
                for (Genome genome : genomes) {
                    log.info("Processing genome {}.", genome);
                    Measurer other = new Measurer(genome, this.usefulRoles);
                    double percent = baseKmers.computePercentSimilarity(other);
                    System.out.format("%s\t%s\t%s\t%s\t%8.2f", baseGenome.getId(),
                            baseGenome.getName(), genome.getId(),
                            genome.getName(), percent);
                    if (this.seedProtein != null) {
                        double seedDist = baseKmers.computeDistance(other, this.seedProtein);
                        System.out.format("\t%8.2f", seedDist);
                    }
                    System.out.println();
                }
            }
        }
    }

}
