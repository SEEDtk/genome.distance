/**
 *
 */
package org.theseed.genome.distance;

import java.util.ArrayList;
import java.util.List;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import org.apache.commons.collections4.map.LRUMap;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.counters.CountMap;
import org.theseed.genome.Genome;
import org.theseed.genome.distance.Measurer.IParms;
import org.theseed.genome.iterator.GenomeSource;
import org.theseed.io.TabbedLineReader;
import org.theseed.p3api.Connection;
import org.theseed.p3api.P3Genome;
import org.theseed.utils.BaseReportProcessor;
import org.theseed.utils.ParseFailureException;

/**
 * This command processes an outlier report to determine the relative performance of SSU- and
 * PheS-based identification of representative genomes.  The outlier report is a tab-delimited
 * file with headers.  Each row contains three different genome IDs-- a source genome ID in
 * "genome_id", a closest genome via PheS similarity in "seed_ref_id", and a closest genome
 * via SSU similarity in "rna_ref_id".  In addition, the PheS-similarity between the source
 * and the SSU-closest is stored in "rna_seed_sim".  The client specifies a minimum threshold
 * for the representative genomes, and if the "rna_seed_sim" value is less than this minimum,
 * the match is considered problematic.  For each problematic match, we compute the full-genome
 * similarity between the source and both representatives.  The resulting percent similarities
 * are output along with an indicator of which is closest.
 *
 * The positional parameters are the name of the genome source for the representative genomes
 * and the representation threshold for the representative set.  The outlier file is
 * presumed to be on the standard input, and the report is written to the standard output.
 *
 * The outlier genomes themselves are read from PATRIC.  We do PROTEINS-level or CONTIGS-level
 * downloads for performance.
 *
 * The command-line options are as follows.
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -i	input file for outlier table (if not STDIN)
 * -t	type of genome source (MASTER, DIR, PATRIC)
 * -R	role definition file for protein comparisons
 * -K	kmer size for contig comparisons
 *
 * --summary		name of a file to contain a summary report; the default is no report
 * --comparison		comparison (PROTEIN or CONTIG); default PROTEIN
 * --cache			cache size for repgen genomes
 *
 * @author Bruce Parrello
 *
 */
public class OutlierCheckProcessor extends BaseReportProcessor implements IParms {

    /**
     * This class contains the useful data on a single genome.
     */
    private class GenomeData {

        /** ID of the test genome */
        private String testGenomeId;
        /** name of the test genome */
        private String testGenomeName;
        /** ID of the PheS representative */
        private String seedGenomeId;
        /** ID of the SSU representative */
        private String rnaGenomeId;

        /**
         * Extract the useful fields from an input line.
         *
         * @param line	input line to process
         */
        protected GenomeData(TabbedLineReader.Line line) {
            this.testGenomeId = line.get(OutlierCheckProcessor.this.testGenomeCol);
            this.testGenomeName = line.get(OutlierCheckProcessor.this.testGenomeNameCol);
            this.seedGenomeId = line.get(OutlierCheckProcessor.this.seedRepCol);
            this.rnaGenomeId = line.get(OutlierCheckProcessor.this.rnaRepCol);
        }

        /**
         * @return the test genome ID
         */
        public String getTestGenomeId() {
            return this.testGenomeId;
        }

        /**
         * @return the test genome name
         */
        public String getTestGenomeName() {
            return this.testGenomeName;
        }

        /**
         * @return the PheS representative ID
         */
        public String getSeedGenomeId() {
            return this.seedGenomeId;
        }

        /**
         * @return the SSU representative ID
         */
        public String getRnaGenomeId() {
            return this.rnaGenomeId;
        }

    }

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(OutlierCheckProcessor.class);
    /** connection to PATRIC */
    private Connection p3;
    /** map of repgen genome IDs to measurers; this is a cache; we can't hold more that 150
     * on most machines */
    private LRUMap<String, Measurer> measureMap;
    /** repgen genome source */
    private GenomeSource repGenomes;
    /** detail level for genomes */
    private P3Genome.Details detailLevel;
    /** input column index for test genome ID */
    private int testGenomeCol;
    /** input column index for test genome name */
    private int testGenomeNameCol;
    /** input column index for seed representative ID */
    private int seedRepCol;
    /** input column index for SSU representative ID */
    private int rnaRepCol;
    /** list of input lines for genomes of interest */
    private List<GenomeData> inputLines;
    /** counts for SSU failures, by genome ID */
    private CountMap<String> rnaFailures;
    /** counts for PheS failures, by genome ID */
    private CountMap<String> seedFailures;
    /** writer for summary report */
    private PrintWriter summStream;
    /** repgen load count */
    private int loadCount;

    // COMMAND-LINE OPTIONS

    /** input file for outlier report (if not STDIN) */
    @Option(name = "--input", aliases = { "-i" }, metaVar = "outlier.tbl", usage = "input file for outlier report (if not STDIN)")
    private File inFile;

    /** type of genome source for reference genomes */
    @Option(name = "--type", aliases = { "-t", "--source" }, usage = "type of repgen genome source")
    private GenomeSource.Type sourceType;

    /** type of comparison to use */
    @Option(name = "--comparison", usage = "comparison algorithm to use")
    private Measurer.Type comparisonType;

    /** size of genome batch */
    @Option(name = "--cache", metaVar = "20", usage = "cache size for loading genomes")
    private int cacheSize;

    /** input role definition file */
    @Option(name = "--roles", aliases = { "-R" }, usage = "role definition file (for method = PROTEIN")
    private File roleFile;

    /** DNA kmer size */
    @Option(name = "--dnaKmer", aliases = { "-K", "--kmer" }, usage = "DNA kmer size (for method = CONTIG")
    private int dnaKmerSize;

    /** summary report file name */
    @Option(name = "--summary", aliases = { "--summ" }, usage = "summary report file name")
    private File summFile;

    /** source for repgen genomes */
    @Argument(index = 0, metaVar = "repgenDir", usage = "genome source for repgen genomes (file or directory)",
            required = true)
    private File repgenDir;

    /** repgen threshold */
    @Argument(index = 1, metaVar = "threshold (int)", usage = "thershold value for the repgen set",
            required = true)
    private int threshold;

    @Override
    protected void setReporterDefaults() {
        this.cacheSize = 100;
        this.sourceType = GenomeSource.Type.DIR;
        this.comparisonType = Measurer.Type.PROTEIN;
    }

    @Override
    protected void validateReporterParms() throws IOException, ParseFailureException {
        // Validate the threshold.
        if (this.threshold < 0)
            throw new ParseFailureException("Threshold must be non-negative.");
        // Validate and build the repgen source
        if (this.cacheSize < 1)
            throw new ParseFailureException("Cache size must be at least 1.");
        if (! this.repgenDir.exists())
            throw new FileNotFoundException("Repgen source " + this.repgenDir + " does not exist.");
        this.repGenomes = this.sourceType.create(this.repgenDir);
        this.measureMap = new LRUMap<String, Measurer>(this.cacheSize);
        this.loadCount = 0;
        // Connect to PATRIC.
        this.p3 = new Connection();
        // Determine the detail level for this comparison type.
        this.detailLevel = this.comparisonType.getLevel();
        // Initialize the comparison parameters.
        this.comparisonType.init(this);
        // Open the input stream.
        TabbedLineReader inStream;
        log.info("Measurement type is {}.", this.comparisonType);
        if (this.inFile == null) {
            log.info("Outlier report will be taken from the standard input.");
            inStream = new TabbedLineReader(System.in);
        } else {
            log.info("Outlier report will be read from {}.", this.inFile);
            inStream = new TabbedLineReader(this.inFile);
        }
        try {
            // Get the columns for the input stream.
            this.testGenomeCol = inStream.findField("genome_id");
            this.testGenomeNameCol = inStream.findField("name");
            this.seedRepCol = inStream.findField("seed_ref_id");
            this.rnaRepCol = inStream.findField("rna_ref_id");
            int threshCol = inStream.findField("rna_seed_sim");
            // Read in the input file and keep the lines of interest.
            this.inputLines = new ArrayList<GenomeData>(2000);
            for (TabbedLineReader.Line line : inStream) {
                int crossSim = (int) line.getDouble(threshCol);
                if (crossSim < this.threshold) {
                    // Here we have an input line of interest.
                    GenomeData lineData = this.new GenomeData(line);
                    this.inputLines.add(lineData);
                }
            }
            log.info("{} outliers of interest found.", this.inputLines.size());
        } finally {
            inStream.close();
        }

        // Finally, insure we can open the output for the summary report.
        if (this.summFile != null) {
            this.summStream = new PrintWriter(this.summFile);
            this.summStream.println("type\tgenome_id\tgenome_name\tfailures");
        }
    }

    @Override
    protected void runReporter(PrintWriter writer) throws Exception {
        try {
            // Initialize the count maps.
            this.seedFailures = new CountMap<String>();
            this.rnaFailures = new CountMap<String>();
            // Write the header line.
            writer.println("genome_id\tname\tseed_rep_prox\trna_rep_prox\tbest\tseed_rep_id\tseed_rep_name\trna_rep_id\trna_rep_name");
            // Start the performance timer.
            long start = System.currentTimeMillis();
            int processCount = 0;
            int seedCount = 0;
            int rnaCount = 0;
            double seedSum = 0.0;
            double rnaSum = 0.0;
            // Loop through the input file.
            for (GenomeData line : this.inputLines) {
                processCount++;
                String testGenomeId = line.getTestGenomeId();
                log.info("Processing genome {}:  {} {}.", processCount, testGenomeId,
                        line.getTestGenomeName());
                Genome testGenome = P3Genome.load(p3, testGenomeId, this.detailLevel);
                // Create the measurer for this genome.
                Measurer measurer = this.comparisonType.create(testGenome);
                // Measure the two representatives against the test genome.
                Measurer seedGenome = this.getMeasurer(line.getSeedGenomeId());
                double seedProximity = measurer.computePercentSimilarity(seedGenome);
                seedSum += seedProximity;
                Measurer rnaGenome = this.getMeasurer(line.getRnaGenomeId());
                double rnaProximity = measurer.computePercentSimilarity(rnaGenome);
                rnaSum += rnaProximity;
                // Figure out if this is a win for PheS or SSU.  The loser gets counted
                // as a failure in the appropriate count map.
                String best = "";
                if (seedProximity > rnaProximity) {
                    best = "PheS";
                    seedCount++;
                    this.rnaFailures.count(rnaGenome.getId());
                } else if (seedProximity < rnaProximity) {
                    best = "SSU";
                    rnaCount++;
                    this.seedFailures.count(seedGenome.getId());
                }
                if (log.isInfoEnabled()) {
                    double rate = (System.currentTimeMillis() - start) / (processCount * 1000.0);
                    log.info("{} genomes compared, {} SSU wins, {} PheS wins, {} loads, {} seconds/genome.",
                            processCount, rnaCount, seedCount, this.loadCount, rate);
                }
                writer.format("%s\t%s\t%6.2f\t%6.2f\t%s\t%s\t%s\t%s\t%s%n", testGenomeId,
                        testGenome.getName(), seedProximity, rnaProximity, best, seedGenome.getId(),
                        seedGenome.getName(), rnaGenome.getId(), rnaGenome.getName());
                // Flush periodically.  This is a slow program, and we need to see progress.
                if (processCount % 10 == 0)
                    writer.flush();
            }
            double seedMean = (seedCount == 0 ? 0 : seedSum / seedCount);
            double rnaMean =  (rnaCount == 0 ? 0 : rnaSum / rnaCount);
            log.info("{} comparisons. {} PheS wins, PheS mean is {}.  {} SSU wins, SSU mean is {}.",
                    processCount, seedCount, seedMean, rnaCount, rnaMean);
            if (this.summStream != null) {
                // Here we need a summary report.
                this.writeSummaryReport("SSU", this.rnaFailures);
                this.writeSummaryReport("PheS", this.seedFailures);
                this.summStream.flush();
            }
        } finally {
            if (this.summStream != null)
                this.summStream.close();
        }
    }

    /**
     * @return a measurer for the specified genome
     *
     * @param genomeId	ID of genome of interest
     * @return
     */
    private Measurer getMeasurer(String genomeId) {
        Measurer retVal = this.measureMap.get(genomeId);
        if (retVal == null) {
            log.info("Loading measurer for genome {}.", genomeId);
            Genome genome = this.repGenomes.getGenome(genomeId);
            retVal = this.comparisonType.create(genome);
            this.measureMap.put(genomeId, retVal);
            this.loadCount++;
        }
        return retVal;
    }

    /**
     * Write out the failure counts for each genome that failed for the specified representation type.
     *
     * @param type		type of failure
     * @param counts	count map for the failures
     */
    private void writeSummaryReport(String type, CountMap<String> counts) {
        for (CountMap<String>.Count counter : counts.sortedCounts()) {
            int count = counter.getCount();
            String genomeId = counter.getKey();
            String genomeName = this.measureMap.get(genomeId).getName();
            this.summStream.format("%s\t%s\t%s\t%d%n", type, genomeId, genomeName, count);
        }
    }

    @Override
    public File getRoleFile() {
        return this.roleFile;
    }

    @Override
    public int getKmerSize() {
        return this.dnaKmerSize;
    }

}
