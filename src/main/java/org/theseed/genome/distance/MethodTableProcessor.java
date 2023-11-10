/**
 *
 */
package org.theseed.genome.distance;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.time.Duration;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import org.apache.commons.lang3.StringUtils;
import org.apache.commons.math3.stat.correlation.KendallsCorrelation;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;
import org.apache.commons.text.TextStringBuilder;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.ParseFailureException;
import org.theseed.counters.CorrelationVariance;
import org.theseed.genome.Genome;
import org.theseed.genome.distance.methods.DistanceMethod;
import org.theseed.genome.distance.methods.GenomePairList;
import org.theseed.genome.distance.methods.Measurer;
import org.theseed.genome.distance.methods.TaxonDistanceMethod;
import org.theseed.genome.iterator.GenomeSource;
import org.theseed.io.TabbedLineReader;
import org.theseed.reports.StringTupleSort;
import org.theseed.utils.BasePipeProcessor;
import org.theseed.utils.StringPair;

/**
 * This command reads an input file describing genome-comparison methods and runs all the methods on pairs
 * of incoming genomes.  The output displays the results of each method on each genome pair.
 *
 * The positional parameters are the name of a file containing the list of methods to test, the name of
 * the role definition file, and the names of one or more genome source directories (or files, in the case of
 * PATRIC sources).
 *
 * The input file is tab-delimited with headers and contains a list of genome ID pairs.  All the genomes must
 * be present in the genome source.
 *
 * The positional parameters are the name of the method list file, the name of a genome source file or directory,
 * and the name of a role definition file for role-parsing.
 *
 * The method list file consists of two columns-- a method type and a parameter string.  The parameter
 * string is free-form and dependent on the method type; generally it consists of space-delimited values.
 *
 * An output file from a previous run can be specified, and if so, the results from the previous run will be used
 * wherever possible.
 *
 * The command-line options are as follows:
 *
 * -h	display command-line usage
 * -v	display more frequent progress messages
 * -i	input file name (if not STDIN)
 * -o	output file name (if not STDOUT)
 * -1	index or name of the input column containing the first genome ID (default "1")
 * -2	index or name of the input column containing the second genome ID (default "2")
 *
 * --source			genome source type (default DIR)
 * --stats			name of an output file for statistics (default "stats.tbl" in the current directory)
 * --previous		if specified, the name of a file containing the results of a previous run to be reused
 *
 * @author Bruce Parrello
 *
 */
public class MethodTableProcessor extends BasePipeProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(MethodTableProcessor.class);
    /** genome sources */
    private List<GenomeSource> genomes;
    /** list of methods to apply */
    private List<DistanceMethod> methods;
    /** list of genome pairs */
    private GenomePairList pairs;
    /** ID of the first genome in the current pair */
    private String genomeId1;
    /** name of the first genome in the current pair */
    private String genomeName1;
    /** taxonomic analysis of the first genome in the current pair */
    private TaxonDistanceMethod.Analysis taxAnalysis1;
    /** taxonomic analyzer */
    private TaxonDistanceMethod taxMethod;
    /** list of distance sets */
    private List<double[]> distanceList;
    /** map of previous results */
    private Map<StringPair, double[]> oldResultMap;

    // COMMAND-LINE OPTIONS

    /** type of genome source */
    @Option(name = "--source", usage = "type of genome source")
    private GenomeSource.Type sourceType;

    /** name of the input column with the first genome ID */
    @Option(name = "--col1", aliases = { "-1", "--c1" }, metaVar = "genome_id",
            usage = "index (1-based) or name of the input column containing the first genome ID")
    private String col1Name;

    /** name of the input column with the second genome ID */
    @Option(name = "--col2", aliases = { "-2", "--c2" }, metaVar = "genome_id",
            usage = "index (1-based) or name of the input column containing the second genome ID")
    private String col2Name;

    /** name of an output file for the correlation statistics */
    @Option(name = "--stats", metaVar = "correlations.tbl", usage = "name of the output file for correlation statistics")
    private File statsFile;

    /** name of the file containing the previous results */
    @Option(name = "--previous", metaVar = "distances.old.tbl",
            usage = "if specified, the name of a file containing previous results that can be re-used")
    private File previousFile;

    /** method list file */
    @Argument(index = 0, metaVar = "methods.tbl", usage = "name of method list file", required = true)
    private File methodFile;

    /** role definition file */
    @Argument(index = 1, metaVar = "roles.in.subsystems", usage = "name of the role definition file", required = true)
    private File roleFile;

    /** genome source directories */
    @Argument(index = 2, metaVar = "inDir1 inDir2 ...", usage = "name of genome source directories (or files)", required = true)
    private List<File> inDirs;

    @Override
    protected void setPipeDefaults() {
        this.sourceType = GenomeSource.Type.DIR;
        this.col1Name = "1";
        this.col2Name = "2";
        this.statsFile = new File("stats.tbl");
        this.previousFile = null;
    }

    @Override
    protected void validatePipeInput(TabbedLineReader inputStream) throws IOException {
        // Here we read the whole input file and form the genome pair list.  First, we need to get the input column IDs.
        int col1idx = inputStream.findField(this.col1Name);
        int col2idx = inputStream.findField(this.col2Name);
        // Create the pair list.
        this.pairs = new GenomePairList();
        // Now read the genome ID pairs.
        log.info("Reading genome pairs.");
        for (var line : inputStream) {
            String g1 = line.get(col1idx);
            String g2 = line.get(col2idx);
            this.pairs.addPair(g1, g2);
        }
        log.info("{} genome pairs to process.", this.pairs.size());
    }

    @Override
    protected void validatePipeParms() throws IOException, ParseFailureException {
        // Set up the role definitions.
        DistanceMethod.loadRoles(this.roleFile);
        // Read in the methods.
        if (! this.methodFile.canRead())
            throw new FileNotFoundException("Method input file " + this.methodFile + " is not found or unreadable.");
        this.methods = new ArrayList<DistanceMethod>();
        int mCount = 0;
        log.info("Reading methods from {}.", this.methodFile);
        try (TabbedLineReader methodStream = new TabbedLineReader(this.methodFile)) {
            for (var line : methodStream) {
                mCount++;
                DistanceMethod method = DistanceMethod.create(line.get(0));
                method.parseParmString(line.get(1));
                log.info("Method {} is {}.", mCount, method);
                this.methods.add(method);
            }
            log.info("{} methods loaded.", mCount);
        }
        // Process the previous-results file.
        if (this.previousFile == null) {
            log.info("No previous results specified.");
            this.oldResultMap = null;
        } else if (! this.previousFile.canRead())
            throw new FileNotFoundException("Previous-results file " + this.previousFile + " is not found or unreadable.");
        else {
            // Here we have previous results to load.  We need to verify that the methods match.
            try (TabbedLineReader prevStream = new TabbedLineReader(this.previousFile)) {
                log.info("Validating previous-results file " + this.previousFile + ".");
                int method0 = prevStream.findField("tax_group") + 1;
                String[] labels = prevStream.getLabels();
                final int nMethods = this.methods.size();
                if (method0 + nMethods != labels.length)
                    throw new IOException("Previous-results file has the wrong number of columns for this method configuration.");
                for (int i = 0; i < this.methods.size(); i++) {
                    if (! labels[i + method0].contentEquals(this.methods.get(i).toString()))
                        throw new IOException("Method " + i + " does not match previous-results file.");
                }
                // Now we know the methods match, so we can read the file.
                this.oldResultMap = new HashMap<StringPair, double[]>(4000);
                log.info("Reading previous results.");
                int id1ColIdx = prevStream.findField("id1");
                int id2Colidx = prevStream.findField("id2");
                for (var line : prevStream) {
                    String id1 = line.get(id1ColIdx);
                    String id2 = line.get(id2Colidx);
                    StringPair pair = new StringPair(id1, id2);
                    double[] distances = new double[nMethods];
                    for (int i = 0; i < nMethods; i++)
                        distances[i] = line.getDouble(method0 + i);
                    this.oldResultMap.put(pair, distances);
                }
                log.info("{} old results read into cache.", this.oldResultMap.size());
            }

        }
        // Connect to the genome sources.
        this.genomes = new ArrayList<GenomeSource>(this.inDirs.size());
        for (File inDir : this.inDirs) {
            if (! inDir.exists())
                throw new FileNotFoundException("Genome source " + inDir + " is not found.");
            this.genomes.add(this.sourceType.create(inDir));
        }
        // Initialize the taxonomic distance method.
        this.taxMethod = new TaxonDistanceMethod();
    }

    @Override
    protected void runPipeline(TabbedLineReader inputStream, PrintWriter writer) throws Exception {
        try {
            // One last validation step:  verify all the genomes are found in the source.
            this.checkGenomes();
            // Prepare the pair list for iteration.
            log.info("Preparing the pair list.");
            this.pairs.prepare();
            // Write the header line.
            writer.println("id1\tname1\tid2\tname2\ttax_group\t" +
                    this.methods.stream().map(x -> x.toString()).collect(Collectors.joining("\t")));
            // Insure we have at least one pair before we continue.
            if (this.pairs.size() > 0) {
                log.info("Initializing method cache.");
                // Get the number of methods.
                final int nMethods = this.methods.size();
                // We will store each pair's distances in here for the correlation statistics later.
                this.distanceList = new ArrayList<double[]>(this.pairs.size());
                // Cache the first pair's first genome.
                this.genomeId1 = this.pairs.get(0).getId1();
                List<Measurer> measurers = this.getMeasurers(genomeId1);
                // Now we loop through the pairs.
                log.info("Processing pairs.");
                int pCount = 0;
                long start = System.currentTimeMillis();
                for (var pair : this.pairs) {
                    double[] distances = new double[nMethods];
                    String id1 = pair.getId1();
                    if (! id1.contentEquals(genomeId1)) {
                        // Here we have a new first genome and we need to re-cache the measurers.
                        this.genomeId1 = id1;
                        measurers = this.getMeasurers(genomeId1);
                    }
                    // Get the second genome.
                    String genomeId2 = pair.getId2();
                    Genome genome = this.getGenome(genomeId2);
                    // Check for a previous result.
                    boolean oldFound = false;
                    if (this.oldResultMap != null)
                        oldFound = this.checkPrevious(this.genomeId1, genomeId2, distances);
                    if (! oldFound) {
                        for (int i = 0; i < nMethods; i++) {
                            DistanceMethod method = this.methods.get(i);
                            // Store the distance.
                            distances[i] = method.getDistance(measurers.get(i), genome);
                        }
                    }
                    // Save the distances.
                    this.distanceList.add(distances);
                    // Compute the taxonomic grouping.
                    final var taxAnalysis2 = this.taxMethod.new Analysis(genome);
                    String taxGroup = this.taxMethod.getGroupingLevel(taxAnalysis1, taxAnalysis2);
                    // Build an output line.
                    TextStringBuilder printLine = new TextStringBuilder(150);
                    printLine.append("%s\t%s\t%s\t%s\t%s", this.genomeId1, this.genomeName1, genomeId2, genome.getName(), taxGroup);
                    for (double distance : distances)
                        printLine.append("\t%8.4f", distance);
                    // Write the line.
                    writer.println(printLine);
                    writer.flush();
                    pCount++;
                    if (log.isInfoEnabled()) {
                        Duration remain = Duration.ofMillis((System.currentTimeMillis() - start) * (this.pairs.size() - pCount) / pCount);
                        log.info("{} pairs processed. {} remaining.", pCount, remain);
                    }
                }
                // Now we compute the statistics.
                try (PrintWriter statWriter = new PrintWriter(this.statsFile)) {
                    this.writeStatistics(statWriter);
                }
            }
        } finally {
            // Here we close the methods to do any needed cleanup (e.g. temporary directories used by BLAST).
            log.info("Cleaning up method storage.");
            for (var method : this.methods)
                method.close();
            this.taxMethod.close();
        }
    }

    /**
     * Check for a previous result in the cache, and use it if one is there.
     *
     * @param id1			first genome ID
     * @param id2			second genome ID
     * @param distances		output array for distances
     *
     * @return TRUE if an old result was found, else FALSE
     */
    private boolean checkPrevious(String id1, String id2, double[] distances) {
        boolean retVal;
        // Create a string pair.
        var pair = new StringPair(id1, id2);
        double[] oldDistances = this.oldResultMap.get(pair);
        if (oldDistances == null)
            retVal = false;
        else {
            // Copy the cached results to the provided result array.
            System.arraycopy(oldDistances, 0, distances, 0, distances.length);
            retVal = true;
        }
        return retVal;
    }

    /**
     * Write the correlation statistics for the method pairs to the specified output.
     *
     * @param statWriter	writer to contain the correlation report
     */
    private void writeStatistics(PrintWriter statWriter) {
        final int n = this.methods.size();
        final int nPairs = this.distanceList.size();
        // These arrays will hold the distances.
        double[] dist1 = new double[nPairs];
        double[] dist2 = new double[nPairs];
        // Create the correlation calculators.
        var pearson = new PearsonsCorrelation();
        var kendall = new KendallsCorrelation();
        var spearman = new SpearmansCorrelation();
        var variance = new CorrelationVariance();
        // Write the header line.
        statWriter.println("method1\tmethod2\tPearson\tKendall\tSpearman\tvariation\tIQR");
        // We generate each correlation line twice, then sort them.
        var sorter = new TreeMap<String[], String>(new StringTupleSort());
        // Loop through every pair of methods.
        for (int i = 0; i < n; i++) {
            String method1 = this.methods.get(i).toString();
            // Form the first method's distance array.
            this.copyDistances(i, dist1);
            for (int j = i+1; j < n; j++) {
                String method2 = this.methods.get(j).toString();
                // Form the second method's distance array.
                this.copyDistances(j, dist2);
                // Compute the correlations.
                var p = pearson.correlation(dist1, dist2);
                var k = kendall.correlation(dist1, dist2);
                var s = spearman.correlation(dist1, dist2);
                var tm = variance.variation(dist1, dist2);
                var iqr = variance.getIQR();
                // Format a line and store it in the sorter for both directions of the comparison.
                String line = String.format("%s\t%s\t%8.4f\t%8.4f\t%8.4f\t%8.4f\t%8.4f", method1, method2, p, k, s, tm, iqr);
                sorter.put(new String[] { method1,  method2 }, line);
                line = String.format("%s\t%s\t%8.4f\t%8.4f\t%8.4f\t%8.4f\t%8.4f", method2, method1, p, k, s, tm, iqr);
                sorter.put(new String[] { method2,  method1 }, line);
            }
        }
        // Unspool the print lines from the sorter.
        sorter.values().stream().forEach(x -> statWriter.println(x));
    }

    /**
     * Copy the distances from the specified method to the specified array.
     *
     * @param i			index of the source method
     * @param dist		target array to build
     */
    private void copyDistances(int i, double[] dist) {
        IntStream.range(0, dist.length).forEach(k -> dist[k] = this.distanceList.get(k)[i]);
    }

    /**
     * This method loads the specified genome and returns a list of measurers for it.
     *
     * @param genomeId	ID of the genome to load
     *
     * @return a list of measurers for the genome, in method order
     */
    private List<Measurer> getMeasurers(String genomeId) {
        List<Measurer> retVal = new ArrayList<Measurer>(this.methods.size());
        Genome genome = this.getGenome(genomeId);
        this.genomeName1 = genome.getName();
        this.taxAnalysis1 = this.taxMethod.new Analysis(genome);
        for (var method : this.methods) {
            var measurer = method.getMeasurer(genome);
            retVal.add(measurer);
        }
        return retVal;
    }

    /**
     * @return the genome with the specified ID
     *
     * @param genomeId	ID of the genome to be harvested from the sources
     */
    private Genome getGenome(String genomeId) {
        Genome retVal = null;
        for (int i = 0; i < this.genomes.size() && retVal == null; i++)
            retVal = this.genomes.get(i).getGenome(genomeId);
        return retVal;
    }

    /**
     * Verify that all the genome IDs in the pair list are in the genome source.
     *
     * @throws IOException
     */
    private void checkGenomes() throws IOException {
        log.info("Validating genome pairs.");
        Set<String> idSet = this.genomes.stream().flatMap(x -> x.getIDs().stream()).collect(Collectors.toSet());
        var errorSet = this.pairs.getIdSet().stream().filter(x -> ! idSet.contains(x)).collect(Collectors.toList());
        if (! errorSet.isEmpty())
            throw new IOException("The following genomes are missing from the sources: " +
                    StringUtils.join(errorSet, ", "));
    }

}
