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
import org.theseed.sequence.GenomeKmers;
import org.theseed.sequence.LSHDiskSeqHash;
import org.theseed.utils.BaseProcessor;

/**
 * This command builds and updates a minhash database of genomes.  The genomes are read from one or more
 * directories and then added to the database.  The database may then be used later to find genomes close
 * to existing genomes.
 *
 * The positional parameters are the name of the database directory followed by the names of the input
 * genome directories.  The command-line options are as follows.
 *
 * -h	display command usage
 * -v	display more detailed progress messages
 * -K 	kmer size for genome signatures (default 21; only used if --create is specified)
 * -w	width of the genome signature (default 2000; only used if --create is specified)
 * -s	number of stages to use (default 20; only used if --create is specified)
 * -b	number of buckets per stage (default 500; only used if --create is specified)
 * -M	optimal number of buckets to keep in memory (default 1000)
 *
 * --create		if specified, it is presumed the database is new and it will be created, erasing all
 * 				current content
 *
 * @author Bruce Parrello
 *
 */
public class BuildProcessor extends BaseProcessor {

    // FIELDS

    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(BuildProcessor.class);
    /** main database */
    private LSHDiskSeqHash genomeHash;

    // COMMAND-LINE OPTIONS

    /** kmer size */
    @Option(name="-K", aliases = { "--kmerSize", "--kmer" }, metaVar = "12", usage = "DNA kmer size")
    private int kmerSize;

    /** width of a genome sketch */
    @Option(name = "-w", aliases = { "--width", "--sketch" }, metaVar = "200", usage = "number of values per genome sketch")
    private int width;

    /** number of stages for hashing */
    @Option(name = "-s", aliases = { "--stages" }, metaVar = "10", usage = "number of hashing stages")
    private int stages;

    /** number of buckets for hashing */
    @Option(name = "-b", aliases = { "--buckets" }, metaVar = "100", usage = "number of hashing buckets")
    private int buckets;

    /** memory cache size */
    @Option(name = "-M", aliases = { "--cache" }, metaVar = "1000", usage = "number of buckets to keep in memory")
    private int cacheLimit;

    /** if specifed, the database will be created */
    @Option(name = "--create", usage = "if specified, the database will be created, erasing the existing directory")
    private boolean createMode;

    /** database directory name */
    @Argument(index = 0, metaVar = "dbDir", usage = "database directory", required = true)
    private File dataDir;

    /** input genome directory names */
    @Argument(index = 1, metaVar = "genomeDir", usage = "input genome directory", multiValued = true)
    private List<File> genomeDirs;


    @Override
    protected void setDefaults() {
        this.width = 2000;
        this.stages = 20;
        this.buckets = 500;
        this.cacheLimit = 1000;
        this.kmerSize = 21;
        this.createMode = false;
        this.genomeHash = null;
    }

    @Override
    protected boolean validateParms() throws IOException {
        if (this.width < 10)
            throw new IllegalArgumentException("Signature width must be 10 or more.");
        if (this.buckets < 1)
            throw new IllegalArgumentException("Buckets per stage must be 1 or more.");
        if (this.stages < 1)
            throw new IllegalArgumentException("Stage count must be 1 or more.");
        if (this.cacheLimit < 1)
            throw new IllegalArgumentException("Memory cache limit must be at least 1.");
        if (this.kmerSize < 1)
            throw new IllegalArgumentException("Invalid kmer size.");
        if (! this.dataDir.isDirectory() && ! this.createMode) {
            throw new FileNotFoundException("Database directory " + this.dataDir + " is not found or invalid.");
        }
        // Validate the genome directories.
        for (File gDir : this.genomeDirs) {
            if (! gDir.isDirectory())
                throw new FileNotFoundException("Genome directory " + gDir + " is not found or invalid.");
        }
        return true;
    }

    @Override
    public void run() {
        try {
            // Set the tuning parameters.
            log.info("Cache size is {} buckets.", this.cacheLimit);
            LSHDiskSeqHash.setCacheLimit(this.cacheLimit);
            if (this.createMode) {
                // Here we have a new database.
                log.info("Creating new database in directory {}.", this.dataDir);
                this.genomeHash = LSHDiskSeqHash.create(this.width, this.stages, this.buckets,
                        this.kmerSize, this.dataDir);
            } else {
                // Here we are loading an existing database.
                log.info("Loading database from directory {}.", this.dataDir);
                this.genomeHash = LSHDiskSeqHash.load(this.dataDir);
            }
            // Set the kmer size.
            GenomeKmers.setKmerSize(this.genomeHash.getKmerSize());
            // Now we loop through the genome directories, loading the genomes.
            for (File gDir : this.genomeDirs) {
                log.info("Processing genome directory {}.", gDir);
                GenomeDirectory genomeDir = new GenomeDirectory(gDir);
                for (Genome genome : genomeDir) {
                    log.info("Loading genome {}.", genome);
                    GenomeKmers kmers = new GenomeKmers(genome);
                    String genomeString = genome.getId() + "\t" + genome.getName();
                    this.genomeHash.add(kmers, genomeString);
                }
            }
            // Save the updated database.
            log.info("Saving genome database.");
            this.genomeHash.save();
        } catch (Exception e) {
            e.printStackTrace();
        } finally {
            // If we have a genome hash, close it.
            if (this.genomeHash != null)
                try {
                    this.genomeHash.close();
                } catch (Exception e) {
                    e.printStackTrace();
                }
        }
    }
}
