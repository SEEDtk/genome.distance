/**
 *
 */
package org.theseed.genome.distance;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.ParseFailureException;
import org.theseed.io.TabbedLineReader;
import org.theseed.sequence.SequenceKmers;
import org.theseed.sequence.hash.Sketch;
import org.theseed.utils.SizeList;

/**
 * This class test different sketch widths for a set of proteins organized into groups.
 * The positional parameters should be the minimum and maximum sketch sizes to try.
 *
 * The error of an observation will be computed as abs(expected - actual) / mean(expected + actual).
 * If the expected distance is 1.0, the observation will be discarded before computing the overall mean
 * error.
 *
 * The input file should be tab-delimited, with group IDs in one column and protein sequences
 * in the other.
 *
 * The command-line options are as follows.
 *
 * -h	display command usage
 * -v	show more detailed progress messages
 * -K	protein kmer size; the default is 8
 * -s	increment for sketch sizes; the default is 10
 * -i	input file containing protein families (default STDIN)
 * -c	index (1-based) or name of the input column containing group IDs
 * -p	index (1-based) or name of the input column containing protein sequences
 * -M	maximum size of a group
 * -e	target mean error level
 *
 * @author Bruce Parrello
 *
 */
public class WidthProcessor extends ProteinKmerReader {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(WidthProcessor.class);
    /** infinity for target size */
    private static final int INVALID_TARGET_SIZE = Integer.MAX_VALUE;
    /** array of sketch sizes to test */
    private int[] sizes;
    /** best sketch size for target error */
    private int targetSize;

    // COMMAND-LINE OPTIONS

    /** sketch size increment */
    @Option(name = "-s", aliases = { "--step", "--incr" }, metaVar = "5", usage = "increment for sketch size search")
    private int stepSize;

    /** maximum group size */
    @Option(name = "-M", aliases = { "--limit", "--maxGroup" }, metaVar = "500", usage = "maximum permissible group size")
    private int maxGroup;

    /** target level for mean error */
    @Option(name = "-e", aliases = { "--error", "--target" }, metaVar = "0.001", usage = "target value for mean error")
    private double targetError;

    /** minimum sketch size to test */
    @Argument(index = 0, metaVar = "50", usage = "starting (minimum) sketch size", required = true)
    private int minSize;

    /** maximum sketch size to test */
    @Argument(index = 1, metaVar = "300", usage = "ending (maximum) sketch size", required = true)
    private int maxSize;

    @Override
    protected void setDefaults() {
        this.stepSize = 10;
        this.initProteinParms();
        this.maxGroup = 1000;
        this.targetError = 0.001;
    }

    @Override
    protected boolean validateParms() throws IOException, ParseFailureException {
        this.validateProteinParms();
        // Validate the positional parameters.
        if (this.minSize > this.maxSize)
            throw new ParseFailureException("Minimum sketch size cannot be larger than maximum.");
        if (this.stepSize <= 0)
            throw new ParseFailureException("Step size must be greater than 0.");
        // Validate the maximum group size.
        if (this.maxGroup < 10)
            throw new ParseFailureException("Maximum group size must be 10 or greater.");
        // Validate the target error level.
        if (this.targetError > 0.1 || this.targetError <= 0.0)
            throw new ParseFailureException("Target error must be > 0 and < 0.1.");
        // Create the size list.
        this.sizes = SizeList.getSizes(this.minSize, this.maxSize, this.stepSize);
        return true;
    }

    /**
     * Process the protein input.
     */
    protected void processProteins() {
        // This will hold the current group ID.
        String groupId = "";
        // The proteins for the current group will accumulate in here.
        List<SequenceKmers> proteins = new ArrayList<SequenceKmers>();
        // This will hold the minimum sketch size to always hit the target.
        this.targetSize = this.minSize;
        // Write the output headers.
        System.out.println("Group\tSize\tPairs\tDwarves\tMean E\tMax E");
        // Loop through the input.
        for (TabbedLineReader.Line line : this.input()) {
            String group = this.getGroupId(line);
            if (! group.contentEquals(groupId) || proteins.size() >= this.maxGroup) {
                // Here we have a new group.  Process the old one if necessary.
                if (proteins.size() > 0)
                    this.ProcessGroup(groupId, proteins);
                // Set up for the new one.
                log.info("Reading group {}.", group);
                groupId = group;
                proteins.clear();
            }
            // Add this protein to the list.
            SequenceKmers prot = this.getProteinKmers(line);
            proteins.add(prot);
        }
        // Process the residual group.  It will be nonempty unless the entire
        // file was empty.
        if (proteins.size() > 0)
            this.ProcessGroup(groupId, proteins);
        // Output the target sketch size.
        if (this.targetSize == WidthProcessor.INVALID_TARGET_SIZE)
            log.warn("Target sketch size is larger than maxmimum.");
        else
            log.info("Target sketch size is {}.", this.targetSize);
    }

    /**
     * Analyze a group.  We test it at each sketch size and output the mean and maximum error rates.
     *
     * @param groupId		name of this group
     * @param sequences		list of the sequences in this group
     */
    private void ProcessGroup(String groupId, List<SequenceKmers> sequences) {
        // First compute the real distances.
        int n = sequences.size();
        log.info("Processing group {} with {} sequences.", groupId, n);
        int pairs = 0;
        double distances[][] = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = i+1; j < n; j++) {
                double distance = sequences.get(i).distance(sequences.get(j));
                if (distance < 1.0) pairs++;
                distances[i][j] = distance;
            }
        }
        if (pairs == 0) {
            log.warn("Group {} has no usable distance pairs.", groupId);
        } else {
            log.info("Group {} has {} usable distance pairs.", groupId, pairs);
            // This will hold the minimum sketch size that meets the target error.
            int minGoodSize = WidthProcessor.INVALID_TARGET_SIZE;
            // Create an array to hold the sketches.
            Sketch[] sketches = new Sketch[n];
            // Now we loop through the sketch sizes.
            for (int size : this.sizes) {
                log.debug("Computing errors for sketch size {}.", size);
                // Fill in the sketches.
                sketches = sequences.stream().map(x -> new Sketch(x.hashSet(size), groupId)).toArray(s -> new Sketch[s]);
                long dwarves = Arrays.stream(sketches).filter(x -> x.getSignature().length < size).count();
                // Loop through the sketches, computing the distance errors.
                double total = 0.0;
                double maxErr = 0.0;
                for (int i = 0; i < n; i++)
                    for (int j = i+1; j < n; j++) {
                        double sketchDist = sketches[i].distance(sketches[j]);
                        double realDist = distances[i][j];
                        if (realDist != sketchDist) {
                            double error = Math.abs(realDist - sketchDist) * 2.0 / (realDist + sketchDist);
                            if (error > maxErr) maxErr = error;
                            total += error;
                        }
                    }
                // Write the results.
                double meanError = total / pairs;
                System.out.format("%s\t%8d\t%8d\t%8d\t%8.4f\t%8.4f%n", groupId, size, pairs, dwarves,
                        meanError, maxErr);
                // Remember it if it's good.
                if (size < minGoodSize && meanError <= this.targetError)
                    minGoodSize = size;
            }
            // Merge in the minimum good size.
            if (minGoodSize > this.targetSize) this.targetSize = minGoodSize;
            if (minGoodSize == WidthProcessor.INVALID_TARGET_SIZE)
                log.warn("{} has no acceptable sketch size in range.", groupId);
            else
                log.info("Minimum acceptable size for {} is {}.", groupId, minGoodSize);
        }
    }

}
