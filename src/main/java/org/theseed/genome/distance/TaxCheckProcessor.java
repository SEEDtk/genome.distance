/**
 *
 */
package org.theseed.genome.distance;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Comparator;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.stream.IntStream;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.excel.utils.Distributor;
import org.theseed.genome.distance.methods.TaxonDistanceMethod;
import org.theseed.io.TabbedLineReader;
import org.theseed.utils.BasePipeProcessor;
import org.theseed.utils.ParseFailureException;

/**
 * This command reads the output from the MethodTableProcessor and computes the mean and standard deviation for each
 * method within each taxonomic grouping.  The taxonomic groupings are in a column labelled "tax_group".  Each
 * subsequent column has the distances for a single method.
 *
 * We will create a map of ranks to statistics objects for each method, then output the results.
 *
 * There are no positional parameters.  The command-line options are as follows.
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -i	input file from MethodTableProcessor (if not STDIN)
 * -o	output file for the report (if not STDOUT)
 *
 * --min	minimum number of data points required to display a result (default 900)
 * --dist	name of an Excel output file for distribution data (default none)
 *
 * @author Bruce Parrello
 *
 */
public class TaxCheckProcessor extends BasePipeProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(TaxCheckProcessor.class);
    /** array of method names */
    private String[] methods;
    /** map of rank names to statistics array; each statistic array parallels "methods" */
    private SortedMap<String, DescriptiveStatistics[]> rankArrays;
    /** distribution tracker */
    private Distributor bucketMap;
    /** index of tax_group column in the input */
    private int taxColIdx;
    /** sorter for the rank names */
    private static final Comparator<String> RANK_SORTER = new TaxonDistanceMethod.RankSort();

    // COMMAND-LINE OPTIONS

    /** minimum number of data points required */
    @Option(name = "--min", aliases = { "-m" }, metaVar = "1000", usage = "minimum number of data points required to display a result")
    private int minPoints;

    /** distribution output file */
    @Option(name = "--dist", metaVar = "buckets.xlsx", usage = "optional output Excel file for distribution data")
    private File distFile;

    @Override
    protected void setPipeDefaults() {
        this.minPoints = 900;
        this.distFile = null;
    }

    @Override
    protected void validatePipeInput(TabbedLineReader inputStream) throws IOException {
        // Find the column containing the taxonomic rank.
        this.taxColIdx = inputStream.findField("tax_group");
        // Get the method list from the remaining columns.
        this.methods = Arrays.copyOfRange(inputStream.getLabels(), this.taxColIdx + 1, inputStream.size());
        // Create the output map.
        this.rankArrays = new TreeMap<String,DescriptiveStatistics[]>(RANK_SORTER);
    }

    @Override
    protected void validatePipeParms() throws IOException, ParseFailureException {
        if (this.minPoints < 1)
            throw new FileNotFoundException("Minimum number of data points must be positive.");
        // Create the distributor.
        this.bucketMap = new Distributor(0.0, 1.0, 50);
    }

    @Override
    protected void runPipeline(TabbedLineReader inputStream, PrintWriter writer) throws Exception {
        // The basic strategy is to read all the inputs into the map and then produce a report on each rank.
        int lCount = 0;
        int method0 = this.taxColIdx + 1;
        for (var line : inputStream) {
            String rank = line.get(this.taxColIdx);
            DescriptiveStatistics[] statMap = this.rankArrays.computeIfAbsent(rank, x -> this.createStats());
            // Now we have an array of summary statistics.  Fill them from the input line.
            for (int i = 0; i < this.methods.length; i++)
                statMap[i].addValue(line.getDouble(method0 + i));
            lCount++;
        }
        log.info("{} data points read from input.", lCount);
        // Now we produce the output.  There is a separate section for each rank, with one line per method and
        // columns for the various statistics.  We start with the header.
        writer.println("rank\tmethod\tmin\tnormal_min\tmean\tnormal_max\tmax\tsdev\tcount");
        // Loop through the methods.
        for (int i = 0; i < this.methods.length; i++) {
            // Loop through the ranks.  They are sorted from biggest (high distance) to smallest (low distance).
            for (var rankEntry : this.rankArrays.entrySet()) {
                String rank = rankEntry.getKey();
                var statArray = rankEntry.getValue();
                String method = this.methods[i];
                var stats = statArray[i];
                // Verify the number of data points.
                int count = (int) stats.getN();
                if (count >= this.minPoints) {
                    // Compute the expected range.
                    double sdev = stats.getStandardDeviation();
                    double mean = stats.getMean();
                    double spread = 2 * sdev;
                    double normalMin = mean - spread;
                    double normalMax = mean + spread;
                    writer.format("%s\t%s\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%d%n",
                            method, rank, stats.getMin(), normalMin, mean, normalMax, stats.getMax(), sdev, stats.getN());
                    // Add this series to the bucket map if we are doing a distribution.
                    if (this.distFile != null) {
                        String seriesName = method + ";" + rank;
                        this.bucketMap.addValues(seriesName, stats.getValues());
                    }
                }
            }
        }
        // Output the distribution if requested.
        if (this.distFile != null)
            this.bucketMap.save(this.distFile);
    }

    /**
     * @return an array of empty summary statistics objects, one per method
     */
    private DescriptiveStatistics[] createStats() {
        DescriptiveStatistics[] retVal = new DescriptiveStatistics[this.methods.length];
        IntStream.range(0, this.methods.length).forEach(i -> retVal[i] = new DescriptiveStatistics());
        return retVal;
    }

}
