/**
 *
 */
package org.theseed.genome.distance;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.theseed.io.TabbedLineReader;
import org.theseed.sequence.ProteinKmers;
import org.theseed.sequence.SequenceKmers;
import org.theseed.sequence.Sketch;
import org.theseed.utils.BaseProcessor;
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
 *
 * @author Bruce Parrello
 *
 */
public class WidthProcessor extends BaseProcessor {

    // FIELDS
    /** input stream */
    private TabbedLineReader inStream;
    /** group ID column */
    private int idIdx;
    /** protein sequence column */
    private int protIdx;
    /** array of sketch sizes to test */
    int[] sizes;

    // COMMAND-LINE OPTIONS

    /** kmer length */
    @Option(name = "-K", aliases = { "--kmer", "--kmerSize" }, metaVar = "12", usage = "protein kmer size")
    private int kmerSize;

    /** sketch size increment */
    @Option(name = "-s", aliases = { "--step", "--incr" }, metaVar = "5", usage = "increment for sketch size search")
    private int stepSize;

    /** input file */
    @Option(name = "-i", aliases = { "--input" }, metaVar = "families.tbl", usage = "input file (if not STDIN)")
    private File inFile;

    /** group ID column name */
    @Option(name = "-c", aliases = { "--col", "--groupCol" }, metaVar = "pgfam_id", usage = "group ID column index (1-based) or name")
    private String idColumn;

    /** protein sequence column name */
    @Option(name = "-p", aliases = { "--prot", "--protCol" }, metaVar = "0", usage = "protein sequence column index (1-based) or name")
    private String protColumn;

    /** minimum sketch size to test */
    @Argument(index = 0, metaVar = "50", usage = "starting (minimum) sketch size", required = true)
    private int minSize;

    /** maximum sketch size to test */
    @Argument(index = 1, metaVar = "300", usage = "ending (maximum) sketch size", required = true)
    private int maxSize;

    @Override
    protected void setDefaults() {
        this.kmerSize = 8;
        this.stepSize = 10;
        this.inFile = null;
        this.idColumn = "1";
        this.protColumn = "aa_sequence";
    }

    @Override
    protected boolean validateParms() throws IOException {
        // Open the input file.
        if (this.inFile == null) {
            log.info("Proteins will be read from standard input.");
            this.inStream = new TabbedLineReader(System.in);
        } else if (! this.inFile.canRead()) {
            throw new FileNotFoundException("Input file " + this.inFile + " is not found or invalid.");
        } else {
            log.info("Proteins will be read from {}.", this.inFile);
            this.inStream = new TabbedLineReader(this.inFile);
        }
        // Find the data columns.
        this.idIdx = this.inStream.findField(this.idColumn);
        this.protIdx = this.inStream.findField(this.protColumn);
        // Validate the positional parameters.
        if (this.minSize > this.maxSize)
            throw new IllegalArgumentException("Minimum sketch size cannot be larger than maximum.");
        if (this.stepSize <= 0)
            throw new IllegalArgumentException("Step size must be greater than 0.");
        // Create the size list.
        this.sizes = SizeList.getSizes(this.minSize, this.maxSize, this.stepSize);
        // Set the kmer size.
        ProteinKmers.setKmerSize(this.kmerSize);
        return true;
    }

    @Override
    public void run() {
        try {
            // This will hold the current group ID.
            String groupId = "";
            // The proteins for the current group will accumulate in here.
            List<SequenceKmers> proteins = new ArrayList<SequenceKmers>();
            // Write the output headers.
            System.out.println("Group\tSize\tPairs\tDwarves\tMean E\tMax E");
            // Loop through the input.
            for (TabbedLineReader.Line line : this.inStream) {
                String group = line.get(this.idIdx);
                if (! group.contentEquals(groupId)) {
                    // Here we have a new group.  Process the old one if necessary.
                    if (proteins.size() > 0)
                        this.ProcessGroup(groupId, proteins);
                    // Set up for the new one.
                    log.info("Reading group {}.", group);
                    groupId = group;
                    proteins.clear();
                }
                // Add this protein to the list.
                SequenceKmers prot = new ProteinKmers(line.get(protIdx));
                proteins.add(prot);
            }
            // Process the residual group.  It will be nonempty unless the entire
            // file was empty.
            if (proteins.size() > 0)
                this.ProcessGroup(groupId, proteins);
            log.info("All done.");
        } catch (Exception e) {
            e.printStackTrace();
        } finally {
            this.inStream.close();
        }
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
                System.out.format("%s\t%8d\t%8d\t%8d\t%8.4f\t%8.4f%n", groupId, size, pairs, dwarves,
                        total / pairs, maxErr);
            }
        }
    }

}
