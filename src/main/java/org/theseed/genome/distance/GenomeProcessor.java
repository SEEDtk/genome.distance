/**
 *
 */
package org.theseed.genome.distance;

import java.io.File;
import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.List;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Genome;
import org.theseed.genome.GenomeDirectory;
import org.theseed.sequence.GenomeKmers;
import org.theseed.utils.BaseProcessor;

/**
 * This class uses kmer distance to compare all the genomes in one or more genome directories using kmers.
 * This is much more memory-intensive than the protein-based comparison. The positional parameters are
 * the name of the directory containing the base GTO files and the name of the genome directory
 * containing the comparison genomes.  Multiple comparison genome directories can be specified.
 *
 * The command-line options are as follows.
 *
 * -v	show more detail on the log
 * -h	display usage information
 * -K	kmer size to use
 *
 * @author Bruce Parrello
 *
 */
public class GenomeProcessor extends BaseProcessor {

    // FIELDS

    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(GenomeProcessor.class);

    // COMMAND-LINE OPTIONS

    /** kmer size */
    @Option(name="-K", aliases = { "--kmerSize", "--kmer" }, metaVar = "12", usage = "DNA kmer size")
    private int kmerSize;

    /** name of the base genome directory */
    @Argument(index=0, metaVar="gtoDir", required=true, usage="base genome GTO directory")
    private File baseDir;

    /** name of the input genome directory */
    @Argument(index=1, metaVar="gtoDir1 gtoDir2 ...", required=true, usage="directory of input GTOs")
    private List<File> genomeDirs;

    @Override
    protected void setDefaults() {
        this.kmerSize = 21;
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
        // Set the kmer size.
        GenomeKmers.setKmerSize(this.kmerSize);
        // We made it this far, we can run the application.
        return true;
    }

    @Override
    public void run() {
        try {
            // Write the output header.
            System.out.println("base_id\tbase_name\tgenome_id\tgenome_name\tdistance");
            // Get the base genome directory.
            GenomeDirectory baseGenomes = new GenomeDirectory(this.baseDir);
            // Create the kmer objects for the base genomes.
            List<GenomeKmers> baseList = scanDirectoryForKmers(baseGenomes);
            // Now loop through the comparison genome directories, computing all the distances.
            for (File compareDir : this.genomeDirs) {
                GenomeDirectory compareGenomes = new GenomeDirectory(compareDir);
                List<GenomeKmers> compareList = this.scanDirectoryForKmers(compareGenomes);
                for (GenomeKmers baseKmers : baseList) {
                    log.info("Processing base genome {} ({}).", baseKmers.getGenomeId(),
                            baseKmers.getGenomeName());
                    for (GenomeKmers compareKmers : compareList) {
                        log.info("Processing compare genome {} ({}).", compareKmers.getGenomeId(),
                                compareKmers.getGenomeName());
                        double dist = baseKmers.distance(compareKmers);
                        System.out.format("%s\t%s\t%s\t%s\t%8.4f%n", baseKmers.getGenomeId(),
                                baseKmers.getGenomeName(), compareKmers.getGenomeId(),
                                compareKmers.getGenomeName(), dist);
                    }
                }
            }
            log.info("All done.");
        } catch (Exception e) {
            e.printStackTrace(System.err);
        }
    }

    /**
     * Read all the genomes from a directory and create a list of kmer objects.
     *
     * @param genomes	directory of genomes to scan
     *
     * @return a list of kmer objects for the genomes
     *
     * @throws NoSuchAlgorithmException
     * @throws UnsupportedEncodingException
     */
    private List<GenomeKmers> scanDirectoryForKmers(GenomeDirectory genomes)
            throws NoSuchAlgorithmException, UnsupportedEncodingException {
        List<GenomeKmers> baseGlist = new ArrayList<GenomeKmers>(genomes.size());
        log.info("Processing genome directory {}.", this.baseDir);
        for (Genome genome : genomes) {
            log.info("Scanning genome {}.", genome);
            GenomeKmers kmers = new GenomeKmers(genome);
            baseGlist.add(kmers);
        }
        return baseGlist;
    }

}
