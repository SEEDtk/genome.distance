/**
 *
 */
package org.theseed.genome.distance;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.SortedSet;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Genome;
import org.theseed.genome.GenomeDirectory;
import org.theseed.sequence.Bucket;
import org.theseed.sequence.GenomeKmers;
import org.theseed.sequence.LSHSeqHash;
import org.theseed.utils.BaseProcessor;

/**
 * This comand generates a sketch database for a directory of genomes.  The sketch database allows
 * estimation of genome distance that is less memory-intensive than the exact comparison used by
 * the "genomes" command.
 *
 * The positional parameters are the name of the query directory and the name of the subject directory.
 * The subject directory will be loaded into a sketch database, and then each query genome will be read
 * in and the closest genomes output out for each.
 *
 * The standard output will be tab-delimited, and contain the ID and name of the query genome, the ID and
 * name of each close subject genome, and the sketch distance.
 *
 * The command-line options are as follows.
 *
 * -h	display command usage
 * -K	kmer size to use
 * -v	display more detailed progress messages
 * -w	width of the sketch database; the default is 2000
 * -s	the number of stages for the locally-sensitive hashing (a higher number is slower but more accurate); the default
 * 		is 15
 * -b	the number of buckets for the locally-sensitive hashing (a lower number is slower but more accurage); the default
 * 		is 100
 * -n	the number of neighbors desired for each query genome; the default is 10
 * -m	maximum acceptable distance for a neighbor genome; the default is 0.9
 *
 * @author Bruce Parrello
 *
 */
public class MashProcessor extends BaseProcessor {

    // FIELDS

    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(MashProcessor.class);
    /** locally-sensitive hash for managing genome sketches */
    private LSHSeqHash subjectHash;
    /** hash of subject genome IDs to names */
    private Map<String, String> genomeNames;

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

    /** number of neighbors to seek */
    @Option(name = "-n", aliases = { "--N", "--neighbors" }, metaVar = "5", usage = "number of close genomes to find")
    private int neighbors;

    /** maximum distance to allow */
    @Option(name = "-m", aliases = { "--max", "--maxDist", "--distance" }, metaVar = "0.75",
            usage = "maximum acceptable distance for a neighboring genome")
    private double maxDist;

    /** directory of query genomes */
    @Argument(index = 0, metaVar = "queryDir", usage = "directory of query genomes", required = true)
    private File queryDir;

    /** directory of subject genomes */
    @Argument(index = 1, metaVar = "subjectDir", usage = "directory of subject genomes", required = true)
    private File subjectDir;


    @Override
    protected void setDefaults() {
        this.width = 2000;
        this.stages = 15;
        this.buckets = 100;
        this.neighbors = 10;
        this.maxDist = 0.9;
        this.kmerSize = 21;
    }

    @Override
    protected boolean validateParms() throws IOException {
        if (! this.queryDir.isDirectory())
            throw new FileNotFoundException("Query directory " + this.queryDir + " not found or invalid.");
        if (! this.subjectDir.isDirectory())
            throw new FileNotFoundException("Subject directory " + this.subjectDir + " not found or invalid.");
        // Create the hash table.
        this.subjectHash = new LSHSeqHash(this.width, this.stages, this.buckets);
        // Set the kmer size.
        GenomeKmers.setKmerSize(this.kmerSize);
        return true;
    }

    @Override
    public void run() {
        try {
            // Initialize some counters.
            int subjectGenomes = 0;
            int queryGenomes = 0;
            int neighborsFound = 0;
            int noneFound = 0;
            // Get the subject directory.
            GenomeDirectory gDir = new GenomeDirectory(this.subjectDir);
            // Create the genome-name hash.
            this.genomeNames = new HashMap<String, String>(gDir.size());
            // Loop through the subject directory, adding genomes to the hash table.
            for (Genome subject : gDir) {
                subjectGenomes++;
                log.info("Processing subject genome #{}: {}.", subjectGenomes, subject);
                GenomeKmers kmers = new GenomeKmers(subject);
                this.genomeNames.put(subject.getId(), subject.getName());
                log.info("Hashing subject genome {}.", subject);
                this.subjectHash.add(kmers, subject.getId());
            }
            log.info("{} subject genomes loaded.", subjectGenomes);
            // We're loaded.  Write the output header.
            System.out.println("query_id\tquery_name\tsubject_id\tsubject_name\tdistance");
            // Now process the query genomes, looking for neighbors.
            gDir = new GenomeDirectory(this.queryDir);
            for (Genome query : gDir) {
                log.info("Searching for neighbors of {}.", query);
                GenomeKmers kmers = new GenomeKmers(query);
                SortedSet<Bucket.Result> results = this.subjectHash.getClosest(kmers, this.neighbors, this.maxDist);
                if (results.size() == 0) {
                    log.warn("No neighbors with distance <= {} found for genome {}.", this.maxDist, query);
                    noneFound++;
                } else {
                    for (Bucket.Result result : results) {
                        // Write out this result.
                        String subjectId = result.getTarget();
                        String subjectName = this.genomeNames.get(subjectId);
                        if (subjectName == null)
                            throw new IllegalStateException("Missing genome " + subjectId + " from subject-genome hash table.");
                        System.out.format("%s\t%s\t%s\t%s\t%8.3f%n", subjectId, subjectName, query.getId(),
                                query.getName(), result.getDistance());
                        neighborsFound++;
                    }
                }
                queryGenomes++;
            }
            log.info("All done. {} genomes processed, {} neighbors found, {} searches failed.",
                    queryGenomes, neighborsFound, noneFound);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

}
