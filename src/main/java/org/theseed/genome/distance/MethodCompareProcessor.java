/**
 *
 */
package org.theseed.genome.distance;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;

import static j2html.TagCreator.*;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.NavigableMap;
import java.util.TreeMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.lang3.StringUtils;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.counters.QualityCountMap;
import org.theseed.genome.Genome;
import org.theseed.genome.distance.methods.Measurer;
import org.theseed.genome.distance.methods.ProtMeasurer;
import org.theseed.genome.iterator.GenomeSource;
import org.theseed.io.TabbedLineReader;
import org.theseed.proteins.RoleMap;
import org.theseed.utils.BasePipeProcessor;
import org.theseed.utils.ParseFailureException;

import j2html.tags.ContainerTag;
import j2html.tags.DomContent;

/**
 * This command compares methods for finding the closest genome.  When the methods
 * disagree, the genomes will be compared using a whole-genome comparison to determine
 * which method was more accurate.
 *
 * The input file is tab-delimited, with headers, and comes in on the standard input.  The
 * first column should contain a genome ID.  Among the subsequent columns are the genome
 * IDs of close genomes.  These are in columns named "XXXXX.closest_genome1", where "XXXXX" is
 * the method name (and may contain periods).  Columns without a header in this format
 * will be ignored.
 *
 * When two methods agree, the agreement is noted.  If they disagree, the chosen close genomes
 * are compared using a whole-genome measurement to determine which is the actual closest
 * genome.
 *
 * The positional parameters are (1) the name of the genome source for the input genomes and (2) the
 * name of the genome source for the representative (close) genomes.  The input genomes occur in the
 * first column and the representative genomes occur in the other columns.
 *
 * The output will be an HTML report on the standard output describing the results.
 *
 * The command-line options are as follows.
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -i	input file (if not STDIN)
 * -o	output report file (if not STDOUT)
 * -t	type of the input genome sources (default DIR)
 * -r	role definition file (default "roles.in.subsystems" in current directory)
 *
 * --repType	type of the representative-genome source (default DIR)
 * --method		method to use for whole-genome comparison (default PROTEIN)
 *
 *
 * @author Bruce Parrello
 *
 */
public class MethodCompareProcessor extends BasePipeProcessor {

    /**
     * This object describes a method failure.
     */
    private static class Failure {

        /** input genome ID */
        private String inGenomeId;
        /** chosen close genome ID */
        private String closeGenomeId;
        /** similarity error */
        private double simError;

        /**
         * Create a failure record.
         *
         * @param inGenome		ID of the input genome with the failure
         * @param closeGenome	ID of the chosen close genome
         * @param error			magnitude of the similarity reduction
         */
        public Failure(String inGenome, String closeGenome, double error) {
            this.inGenomeId = inGenome;
            this.closeGenomeId = closeGenome;
            this.simError = error;
        }

        /**
         * @return the similarity reduction
         */
        public double getSimError() {
            return this.simError;
        }

        /**
         * @return a table header for a failure table
         */
        public static DomContent getHeader() {
            return tr(th("Input Genome"), th("Chosen Rep"), th("Error"));
        }

        /**
         * @return a table row for this failure
         */
        public DomContent getTableRow() {
            return tr(td(this.inGenomeId), td(this.closeGenomeId), td(String.format("%6.2f", this.simError)));
        }

    }

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(MethodCompareProcessor.class);
    /** input genome source */
    private GenomeSource inputGenomes;
    /** representative genome source */
    private GenomeSource repGenomes;
   /** winner and loser counts for each method */
    private QualityCountMap<String> counters;
    /** map of method names to input column indices */
    private Map<String, Integer> colMap;
    /** map of methods to failures */
    private Map<String, List<Failure>> failureMap;
    /** match pattern for closest-genome columns */
    private static final Pattern COLUMN_HEAD = Pattern.compile("(\\S+)\\.closest_genome1");
    /** style sheet for web page */
    private static final String CSS_LINK = "https://core.theseed.org/SEEDtk/css/erdb.css";

    // COMMAND-LINE OPTIONS

    /** type of input genome source */
    @Option(name = "--type", aliases = { "-t" }, usage = "source type for input genomes")
    private GenomeSource.Type inputType;

    /** type of representative-genome source */
    @Option(name = "--repType", usage = "source type for representative genomes")
    private GenomeSource.Type repType;

    /** type of whole-genome comparison method */
    @Option(name = "--method", usage = "whole-genome comparison method")
    private Measurer.Type mMethod;

    /** role definition map file */
    @Option(name = "--roles", aliases = { "-r" }, metaVar = "roles.in.subsystems",
            usage = "name of role definition file (for PROTEIN comparison method)")
    private File roleFile;

    /** file or directory containing input genomes */
    @Argument(index = 0, metaVar = "inDir", usage = "source directory or file for input genomes", required = true)
    private File inDir;

    /** file or directory containing representative genomes */
    @Argument(index = 1, metaVar = "repDir", usage = "source directory or file for representative genomes",
            required = true)
    private File repDir;

    @Override
    protected void setPipeDefaults() {
        this.inputType = GenomeSource.Type.DIR;
        this.repType = GenomeSource.Type.DIR;
        this.mMethod = Measurer.Type.PROTEIN;
        this.roleFile = new File(System.getProperty("cur.dir"), "roles.in.subsystems");
    }

    @Override
    protected void validatePipeInput(TabbedLineReader inputStream) throws IOException {
        // Here we determine which methods are in which columns.
        this.colMap = new TreeMap<String, Integer>();
        String[] labels = inputStream.getLabels();
        for (int i = 1; i < labels.length; i++) {
            Matcher m = COLUMN_HEAD.matcher(labels[i]);
            if (m.matches()) {
                // Here we have the column for a method.
                this.colMap.put(m.group(1), i);
            }
        }
        log.info("{} methods found: {}.", this.colMap.size(), StringUtils.join(this.colMap.keySet(), ", "));
    }

    @Override
    protected void validatePipeParms() throws IOException, ParseFailureException {
        // Verify that the input directories are valid.
        if (! this.inDir.exists())
            throw new FileNotFoundException("Input genome source " + this.inDir + " is not found.");
        if (! this.repDir.exists())
            throw new FileNotFoundException("Representative genome source " + this.inDir + " is not found.");
        // Check for the role file.
        if (! this.roleFile.canRead()) {
            if (this.mMethod == Measurer.Type.PROTEIN)
                throw new FileNotFoundException("A valid role definition file is required for PROTEIN comparison.");
        } else {
            log.info("Loading role definitions from {}.", this.roleFile);
            ProtMeasurer.setRoleMap(RoleMap.load(this.roleFile));
        }
        // Now create the genome sources.
        log.info("Connecting to input genome source.");
        this.inputGenomes = this.inputType.create(this.inDir);
        log.info("Connecting to representative-genome source.");
        this.repGenomes = this.repType.create(this.repDir);
        // Announce the measurer.
        log.info("Whole-genome measurement style is {}.", this.mMethod);
    }

    @Override
    protected void runPipeline(TabbedLineReader inputStream, PrintWriter writer) throws Exception {
        // Create the failure map.
        this.failureMap = new TreeMap<String, List<Failure>>();
        // Create the quality counter.
        this.counters = new QualityCountMap<String>();
        // Now we loop through the input genomes.
        int count = 0;
        for (TabbedLineReader.Line line : inputStream) {
            // For the current input genome, this will map each method to its chosen genome.
            Map<String, String> methodToClosest = new TreeMap<String, String>();
            // For the current input genome, this will map each representative genome to its similarity score.
            NavigableMap<String, Double> genomeToSim = new TreeMap<String, Double>();
            // Get the input genome ID.
            String inGenomeId = line.get(0);
            count++;
            log.info("Processing genome {}: {}.", count, inGenomeId);
            // Get the closest genomes for the methods.
            for (Map.Entry<String, Integer> methodEntry : this.colMap.entrySet()) {
                String method = methodEntry.getKey();
                String closeGenome = line.get(methodEntry.getValue());
                if (! StringUtils.isBlank(closeGenome)) {
                    methodToClosest.put(method, closeGenome);
                    genomeToSim.put(closeGenome, 0.0);
                }
            }
            // Now we need to compute the best genome of the candidates.
            String bestGenome = "";
            double bestSim = 0.0;
            switch (genomeToSim.size()) {
            case 0:
                log.warn("No close genomes found for {}.", inGenomeId);
                break;
            case 1:
                bestGenome = genomeToSim.firstKey();
                break;
            default:
                // Here we need to compare the genomes.
                log.info("Resolving conflict for {}.", inGenomeId);
                Genome inGenome = this.inputGenomes.getGenome(inGenomeId);
                if (inGenome == null)
                    throw new IOException("Input genome " + inGenomeId + " not found in " + this.inDir + ".");
                Measurer measure = this.mMethod.create(inGenome);
                for (Map.Entry<String, Double> closeEntry : genomeToSim.entrySet()) {
                    String repGenomeId = closeEntry.getKey();
                    Genome repGenome = this.repGenomes.getGenome(repGenomeId);
                    if (repGenome == null)
                        throw new IOException("Representative genome " + repGenomeId + " not found in " + this.repDir + ".");
                    double sim = measure.computePercentSimilarity(repGenome);
                    closeEntry.setValue(sim);
                    if (sim > bestSim) {
                        bestGenome = repGenomeId;
                        bestSim = sim;
                    }
                }
            }
            // Now compute the success of each method.
            for (Map.Entry<String, String> methodEntry : methodToClosest.entrySet()) {
                String method = methodEntry.getKey();
                String closeGenome = methodEntry.getValue();
                if (closeGenome.equals(bestGenome))
                    this.counters.setGood(method);
                else {
                    this.counters.setBad(method);
                    if (log.isInfoEnabled()) {
                        double error = bestSim - genomeToSim.get(closeGenome);

                        log.info("Best genome for {} was {}, but method {} chose {}, which is {} less similar.",
                                inGenomeId, bestGenome, method, closeGenome, error);
                        // Add the failure to the failure map.
                        Failure failure = new Failure(inGenomeId, bestGenome, error);
                        List<Failure> failList = this.failureMap.computeIfAbsent(method,
                                x -> new ArrayList<Failure>(this.colMap.size()));
                        failList.add(failure);
                    }
                }
            }
        }
        // Now we build our output.  First, we create a div containing the failure tables for each method,
        // accumulating the mean error along the way.
        ContainerTag methodTables = div();
        var errorMap = new TreeMap<String, Double>();
        for (Map.Entry<String, List<Failure>> failureData : this.failureMap.entrySet()) {
            String method = failureData.getKey();
            List<Failure> failures = failureData.getValue();
            double errorSum = 0.0;
            int errorCount = 0;
            methodTables.with(h3(method));
            ContainerTag failTable = table().with(Failure.getHeader());
            for (Failure failure : failures) {
                failTable.with(failure.getTableRow());
                errorSum += failure.getSimError();
                errorCount++;
            }
            // Record the mean error.
            errorMap.put(method, errorSum / errorCount);
            // Put the method table in the div.
            methodTables.with(failTable);
        }
        // Now we create the method summary table.
        ContainerTag summTable = table().with(tr(th("Method"), th("Success"), th("Failure"), th("Percent"), th("Mean Error")));
        for (String method : this.counters.bestKeys()) {
            int good = this.counters.good(method);
            int bad = this.counters.bad(method);
            String pct = String.format("%6.2f", this.counters.fractionGood(method) * 100.0);
            String error = String.format("%6.3f", errorMap.getOrDefault(method, 0.0));
            summTable.with(tr(td(method), td(Integer.toString(good)).withClass("num"),
                    td(Integer.toString(bad)).withClass("num"), td(pct).withClass("num"),
                    td(error).withClass("num")));
        }
        // Assemble the page body.
        ContainerTag body = body().with(h1("Close-Genome Method Comparison"), h2("Summary"),
                summTable);
        if (errorMap.size() > 0)
            body.with(h2("Method Failures"), methodTables);
        // Now build and output the page.
        ContainerTag page = html().with(head().with(title("Close-Genome Methods")),
                link().withRel("styleSheet").withHref(CSS_LINK)).with(body);
        writer.println(page.renderFormatted());
    }

}
