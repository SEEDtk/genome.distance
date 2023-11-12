/**
 *
 */
package org.theseed.genome.distance;

import java.io.IOException;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.ParseFailureException;
import org.theseed.counters.CountMap;
import org.theseed.sequence.ProteinKmers;

/**
 * This subcommand counts how many times each kmer occurs in the specified set of proteins.  It then outputs the number of
 * kmers found, the mean count value, the maximum count value, and the standard deviation to the log.
 *
 * A maximum number of kmers is specified to insure we produce output before we would overflow the hash table of counts.
 *
 * The command-line options are as follows:
 *
 * -h	display command usage
 * -v	show more detailed progress messages
 * -K	protein kmer size; the default is 8
 * -i	input file containing protein families (default STDIN)
 * -c	index (1-based) or name of the input column containing group IDs
 * -p	index (1-based) or name of the input column containing protein sequences
 *
 * --max	maximum number of kmers to count (default 1 billion)
 *
 * @author Bruce Parrello
 *
 */
public class KmerCountProcessor extends ProteinKmerReader {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(KmerCountProcessor.class);
    /** count map for kmers */
    private CountMap<String> kmerCounts;

    // COMMAND-LINE OPTIONS

    /** maximum number of kmers allowed */
    @Option(name = "--max", metaVar = "900000000", usage = "maximum number of kmers to count")
    private int maxKmers;

    @Override
    protected void setDefaults() {
        this.maxKmers = 1000000000;
        this.initProteinParms();
    }

    @Override
    protected boolean validateParms() throws IOException, ParseFailureException {
        this.validateProteinParms();
        // Insure the kmer max is valid.
        if (this.maxKmers < 10)
            throw new ParseFailureException("Max kmers must be 10 or more.");
        return true;
    }

    @Override
    protected void processProteins() throws IOException {
        // Initialize the count map.
        this.kmerCounts = new CountMap<String>();
        // Loop through the input file.
        int protCount = 0;
        int skipCount = 0;
        long lastMsg = System.currentTimeMillis();
        for (var line : this.input()) {
            protCount++;
            // Get the kmers for this protein.
            ProteinKmers kmers = this.getProteinKmers(line);
            for (String kmer : kmers) {
                if (this.kmerCounts.size() < this.maxKmers)
                    this.kmerCounts.count(kmer);
                else {
                    var count = this.kmerCounts.findCounter(kmer);
                    if (count == null)
                        skipCount++;
                    else
                        count.increment();
                }
            }
            if (log.isInfoEnabled() && System.currentTimeMillis() - lastMsg >= 10000) {
                log.info("{} proteins processed.  {} kmers skipped.  {} in table.", protCount, skipCount, this.kmerCounts.size());
                lastMsg = System.currentTimeMillis();
            }
        }
        log.info("Processing kmer counts.  {} proteins processed, {} kmers skipped, {} kmers found.", protCount, skipCount,
                this.kmerCounts.size());
        SummaryStatistics stats = new SummaryStatistics();
        for (var count : this.kmerCounts.counts())
            stats.addValue(count.getCount());
        log.info("Maximum kmer count is {}, mean is {}, standard deviation is {}.", stats.getMax(), stats.getMean(), stats.getStandardDeviation());
    }


}
