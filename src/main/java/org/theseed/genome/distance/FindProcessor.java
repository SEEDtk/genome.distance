/**
 *
 */
package org.theseed.genome.distance;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.List;
import java.util.SortedSet;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.BaseProcessor;
import org.theseed.basic.ParseFailureException;
import org.theseed.genome.Genome;
import org.theseed.genome.GenomeDirectory;
import org.theseed.sequence.GenomeKmers;
import org.theseed.sequence.hash.Bucket;
import org.theseed.sequence.hash.LSHDiskSeqHash;

/**
 * This command reads all the genomes from a directory and returns the close genomes found in a minHash database
 * created using the BuildProcessor.  The positional parameters are the name of the database directory and the
 * names of the input genome directories.  The results will be written to the standard output.
 *
 * The command-line options are as follows.
 *
 * -h	display command usage
 * -v	display more detailed progress messages
 * -M	optimal number of buckets to keep in memory (default 1000)
 * -n	the number of neighbors desired for each query genome; the default is 10
 * -m	maximum acceptable distance for a neighbor genome; the default is 0.9
 *
 * @author Bruce Parrello
 *
 */
public class FindProcessor extends BaseProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(FindProcessor.class);

    // COMMAND-LINE OPTIONS

    /** memory cache size */
    @Option(name = "-M", aliases = { "--cache" }, metaVar = "1000", usage = "number of buckets to keep in memory")
    private int cacheLimit;

    /** number of neighbors to seek */
    @Option(name = "-n", aliases = { "--N", "--neighbors" }, metaVar = "5", usage = "number of close genomes to find")
    private int neighbors;

    /** maximum distance to allow */
    @Option(name = "-m", aliases = { "--max", "--maxDist", "--distance" }, metaVar = "0.75",
            usage = "maximum acceptable distance for a neighboring genome")
    private double maxDist;

    /** database directory name */
    @Argument(index = 0, metaVar = "dbDir", usage = "database directory", required = true)
    private File dataDir;

    /** input genome directory names */
    @Argument(index = 1, metaVar = "genomeDir", usage = "input genome directory", multiValued = true)
    private List<File> genomeDirs;

    @Override
    protected void setDefaults() {
        this.cacheLimit = 1000;
        this.neighbors = 10;
        this.maxDist = 0.9;
    }

    @Override
    protected boolean validateParms() throws IOException, ParseFailureException {
        if (this.cacheLimit < 1)
            throw new ParseFailureException("Memory cache size must be greater than 0.");
        if (! this.dataDir.isDirectory())
            throw new FileNotFoundException("Genome database directory " + this.dataDir + " is not found or invalid.");
        for (File gDir : this.genomeDirs)
            if (! gDir.isDirectory())
                throw new FileNotFoundException("Genome input directory " + gDir + " is not found or invalid.");
        return true;
    }

    @Override
    public void runCommand() throws Exception {
        // Set the cache limit.
        log.info("Cache limit is {} buckets.", this.cacheLimit);
        LSHDiskSeqHash.setCacheLimit(this.cacheLimit);
        // Load in the genome database.
        try (LSHDiskSeqHash genomeDb = LSHDiskSeqHash.load(this.dataDir)) {
            // Set the genome kmer size.
            log.info("Genome kmer size is {}.", genomeDb.getKmerSize());
            GenomeKmers.setKmerSize(genomeDb.getKmerSize());
            // Write the output header.  Note the neighbor ID and name will come back as one field
            // from the genome database, with an internal delimiting tab.
            System.out.println("genome_id\tgenome_name\tneighbor_id\tneighbor_name\tdistance");
            // Loop through the input directories.
            int foundCount = 0;
            int failCount = 0;
            for (File gDir : this.genomeDirs) {
                log.info("Processing input directory {}.", gDir);
                GenomeDirectory genomes = new GenomeDirectory(gDir);
                for (Genome genome : genomes) {
                    log.info("Processing genome {}.", genome);
                    GenomeKmers kmers = new GenomeKmers(genome);
                    SortedSet<Bucket.Result> neighbors = genomeDb.getClosest(kmers, this.neighbors, this.maxDist);
                    if (neighbors.size() == 0) failCount++;
                    for (Bucket.Result neighbor : neighbors) {
                        System.out.format("%s\t%s\t%s\t%8.3f%n", genome.getId(), genome.getName(),
                                neighbor.getTarget(), neighbor.getDistance());
                        foundCount++;
                    }
                }
            }
            log.info("All done. {} neighbors found. {} failures.", foundCount, failCount);
        }

    }

}
