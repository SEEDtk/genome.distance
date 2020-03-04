/**
 *
 */
package org.theseed.genome.distance;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.theseed.genome.Genome;
import org.theseed.genome.GenomeDirectory;
import org.theseed.proteins.RoleMap;
import org.theseed.utils.BaseProcessor;

/**
 * Compare base genomes to all the genomes in a directory.  The output contains
 * the percent similarity for each input genome by genome ID.  The positional parameters are
 * the name of the directory containing the base GTO files and the name of the genome directory
 * containing the comparison genomes.  Multiple comparison genome directories can be specified.
 *
 * @author Bruce Parrello
 *
 */
public class DistanceProcessor extends BaseProcessor {

    // FIELDS
    /** map of useful roles */
    RoleMap usefulRoles;
    /** current input genome directory */
    GenomeDirectory gtoDir;

    // COMMAND LINE

    /** name of the useful-role file */
    @Option(name="-R", aliases={"--roles"}, metaVar="roleFile", required=true,
            usage="file of useful roles")
    private File roleFile;

    /** name of the base genome directory */
    @Argument(index=0, metaVar="gtoDir", required=true, usage="base genome GTO directory")
    private File baseDir;

    /** name of the input genome directory */
    @Argument(index=1, metaVar="gtoDir1 gtoDir2 ...", required=true, usage="directory of input GTOs")
    private List<File> genomeDirs;

    @Override
    protected void setDefaults() {
    }

    @Override
    protected boolean validateParms() throws IOException {
        // Get the role file.
        this.usefulRoles = RoleMap.load(this.roleFile);
        // Verify the base genome directory.
        if (! this.baseDir.isDirectory())
            throw new IOException("Directory " + this.baseDir + " not found or invalid.");
        // Verify the genome directories.
        for (File genomeDir : this.genomeDirs) {
            if (! genomeDir.isDirectory()) {
                throw new IOException("Directory " + genomeDir + " not found or invalid.");
            }
        }
        // We made it this far, we can run the application.
        return true;
    }

    public void run() {
        try {
            // Write the output header.
            System.out.println("base_id\tbase_name\tgenome_id\tgenome_name\tsimilarity");
            // Create the main measurement object.
            GenomeDirectory baseGenomes = new GenomeDirectory(this.baseDir);
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
                        double percent = baseKmers.computePercentSimilarity(genome);
                        System.out.format("%s\t%s\t%s\t%s\t%8.2f%n", baseGenome.getId(),
                                baseGenome.getName(), genome.getId(),
                                genome.getName(), percent);
                    }
                }
            }
        } catch (IOException e) {
            System.err.println(e.getMessage());
        }
    }

}
