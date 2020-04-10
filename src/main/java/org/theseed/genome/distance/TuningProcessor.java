/**
 *
 */
package org.theseed.genome.distance;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Set;

import org.slf4j.LoggerFactory;
import org.theseed.sequence.Bucket;
import org.theseed.sequence.LSHSeqHash;
import org.theseed.sequence.Sketch;
import org.theseed.utils.BaseProcessor;
import org.theseed.utils.SizeList;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;

/**
 * This command tries different stage and bucket sizes for a protein LSHSeqHash and outputs a qaulity measure.
 * The quality measure is based on a target distance.  During the initial pass over the buckets, we compute the
 * number of sketches within a certain distance of each other sketch and save this in a count map.  The sketch's
 * group ID is replaced by a unique name string so we can map the sketch uniquely.  The quality is the percent
 * of close sketches found for each input sketch.
 *
 * The positional parameters are the name of the bucket file, the minimum stage count, and the maximum stage
 * count.  The command-line options are as follows.
 *
 * -h	display usage
 * -v	display more progress messages on the log
 * -b	number of buckets per stage (default 300)
 * -s	increment between stage counts (default 5)
 * -w	sketch width (default 360)
 * -t	target distance (default 0.7)
 *
 *
 * @author Bruce Parrello
 *
 */
public class TuningProcessor extends BaseProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(TuningProcessor.class);
    /** list of stage sizes to check */
    private int[] stageSizes;

    // COMMAND-LINE OPTIONS

    /** number of buckets per stage */
    @Option(name = "-b", aliases = { "--buckets" }, metaVar = "200", usage = "number of buckets per stage")
    private int bucketCount;

    /** step size between stage counts */
    @Option(name = "-s", aliases = { "--step", "--incr" }, metaVar = "5", usage = "increment for stage count search")
    private int stepSize;

    /** width of a protein sketch */
    @Option(name = "-w", aliases = { "--width", "--sketch" }, metaVar = "200", usage = "number of values per protein sketch")
    private int width;

    /** target distance for a close protein */
    @Option(name = "-t", aliases = { "--target", "--minDist" }, metaVar = "0.50", usage = "target sketch distance")
    private double target;

    /** name of input file */
    @Argument(index = 0, metaVar = "sketchesIn.ser", usage = "input file containing protein sketches", required = true)
    private File inFile;

    /** minimum stage count */
    @Argument(index = 1, metaVar = "10", usage = "starting (minimum) stage count", required = true)
    private int minStageCount;

    /** maximum stage count */
    @Argument(index = 2, metaVar = "100", usage = "ending (maximum) stage count", required = true)
    private int maxStageCount;

    @Override
    protected void setDefaults() {
        this.bucketCount = 300;
        this.stepSize = 5;
        this.width = 360;
        this.target = 0.7;
    }

    @Override
    protected boolean validateParms() throws IOException {
        if (! this.inFile.canRead())
            throw new FileNotFoundException("Input file " + this.inFile + " not found or unreadable.");
        if (this.minStageCount < 1)
            throw new IllegalArgumentException("Minimum stage count must be at least 1.");
        if (this.maxStageCount < this.minStageCount)
            throw new IllegalArgumentException("Maximum stage count must be no less than minimum.");
        if (this.stepSize < 1)
            throw new IllegalArgumentException("Step size must be at least 1.");
        if (this.bucketCount < 10)
            throw new IllegalArgumentException("Bucket count must be at least 10.");
        if (this.target <= 0.0 || this.target >= 1.0)
            throw new IllegalArgumentException("Target distance must be between 0 and 1 (exclusive).");
        // Create the size list.
        this.stageSizes = SizeList.getSizes(this.minStageCount, this.maxStageCount, this.stepSize);
        return true;
    }

    @Override
    public void run() {
        try {
            // Write the output header.
            System.out.println("Stages\tFound\tFailed\tQuality");
            // Read in all the sketches.
            log.info("Reading sketches from {}.", this.inFile);
            Bucket testSketches = Bucket.load(this.inFile);
            int n = testSketches.size();
            log.info("{} proteins found in file.", n);
            // Give each sketch a unique name.
            int idx = 1;
            for (Sketch sketch : testSketches)
                sketch.setName(String.format("p%d", idx++));
            // This bucket will hold the sketches that have close neighbors.
            Bucket goodSketches = new Bucket();
            // Now track the number of close sketches.  Every sketch is trivially close to itself,
            // so we count that as well.
            int totalPairs = 0;
            for (int i = 0; i < n; i++) {
                Sketch sketch1 = testSketches.get(i);
                log.debug("Counting pairs for {}.", sketch1.getName());
                int expected = 0;
                for (Sketch sketch2 : testSketches.after(i)) {
                    double distance = sketch1.distance(sketch2);
                    if (distance < this.target)
                        expected++;
                }
                if (expected > 0) {
                    totalPairs += expected;
                    goodSketches.add(sketch1);
                }
            }
            log.info("{} close pairs found in protein list. {} sequences have neighbors.", totalPairs, goodSketches.size());
            // The expected total is twice the number of pairs (once each direction).
            totalPairs += totalPairs;
            // Loop through the sizes, testing each one.  Note the width doesn't matter, since we are going
            // sketch-to-sketch.
            for (int stageSize : this.stageSizes) {
                log.info("Testing {} stages.", stageSize);
                LSHSeqHash hash = new LSHSeqHash(200, stageSize, this.bucketCount);
                for (Sketch sketch : testSketches)
                    hash.add(sketch);
                log.info("Hash loaded with {} proteins.", testSketches.size());
                // Run through each sketch, finding close sketches.  We subtract one for
                // the sketch that finds itself.
                int found = 0;
                int failed = 0;
                for (Sketch sketch : goodSketches) {
                    Set<Bucket.Result> results = hash.getClose(sketch, this.target);
                    found += results.size() - 1;
                    if (results.size() <= 1)
                        failed++;
                }
                System.out.format("%8d\t%8d\t%8d\t%8.4f%n", stageSize, found, failed,
                        ((double) found) / totalPairs);
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

}
