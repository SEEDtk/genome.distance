/**
 *
 */
package org.theseed.genome.signatures;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.reports.SignatureReporter;
import org.theseed.genome.Genome;
import org.theseed.genome.iterator.GenomeSource;
import org.theseed.basic.BaseReportProcessor;
import org.theseed.basic.ParseFailureException;
import org.theseed.counters.CountMap;

/**
 * This command looks for protein signatures in genome groups.  Each protein is identified by a protein class.  It is possible
 * for a protein to be in more than one class, though this is rare.  (Here we are thinking of the case where the class is a
 * role.)  The incoming genomes are divided into an IN group and an OUT group.  The protein signature consists of the
 * classes that are common in the IN group and uncommon in the OUT group.
 *
 * The positional parameters are the names of the genome sources for the IN group and the OUT group.  The command-line options
 * are as follows.
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -o	output file for the report (if not STDIN)
 *
 * --format		output format for the report (default COUNTS)
 * --t1			source type for the IN group (default DIR)
 * --t2			source type for the OUT group (default DIR)
 * --class		type of protein classification to use (default PGFAM)
 * --min1		minimum fraction of occurrences for a classification in the IN group (default 0.8)
 * --max2		maximum fraction of occurrences for a classification in the OUT group (default 0.2)
 * --roles		name of role definition file (only required for class ROLE)
 * --both		show signatures for both groups
 *
 * @author Bruce Parrello
 *
 */
public class SignatureProcessor extends BaseReportProcessor implements SignatureClass.IParms, SignatureReporter.IParms {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(SignatureProcessor.class);
    /** signature classifier */
    private SignatureClass classifier;
    /** reporting object */
    private SignatureReporter reporter;
    /** IN-group genoems */
    private GenomeSource genomes1;
    /** OUT-group genomes */
    private GenomeSource genomes2;

    // COMMAND-LINE OPTIONS

    /** format for the report */
    @Option(name = "--format", usage = "output report format")
    private SignatureReporter.Type outFormat;

    /** source type for the IN group */
    @Option(name = "--t1", usage = "source type for first (IN) genome group")
    private GenomeSource.Type sourceType1;

    /** source type for the OUT group */
    @Option(name = "--t2", usage = "source type for second (OUT) genome group")
    private GenomeSource.Type sourceType2;

    /** classification type */
    @Option(name = "--class", usage = "classification type")
    private SignatureClass.Type classType;

    /** minimum fraction of occurrences for a classification in the IN group */
    @Option(name = "--min1", usage = "minimum fraction of IN-group genomes that must contain a classification")
    private double min1;

    /** maximum fraction of occurrences for a classification in the OUT group */
    @Option(name = "--max1", usage = "maximum fraction of OUT-group genomes that can contain a classification")
    private double max2;

    /** role definition file */
    @Option(name = "--roles", usage = "role definition file (for class type ROLE")
    private File roleFile;

    /** compute signatures in both directions */
    @Option(name = "--both", usage = "if specified, signatures will be computed in both directions")
    private boolean bothFlag;

    /** input genome source 1 (IN) */
    @Argument(index = 0, metaVar = "genomeDirIN", usage = "source for first (IN) group of genomes")
    private File genomeDir1;

    /** input genome source 2 (OUT) */
    @Argument(index = 1, metaVar = "genomeDirOUT", usage = "source for second (OUT) group of genoems")
    private File genomeDir2;

    @Override
    protected void setReporterDefaults() {
        this.outFormat = SignatureReporter.Type.COUNTS;
        this.sourceType1 = GenomeSource.Type.DIR;
        this.sourceType2 = GenomeSource.Type.DIR;
        this.classType = SignatureClass.Type.PGFAM;
        this.min1 = 0.80;
        this.max2 = 0.20;
        this.roleFile = null;
        this.bothFlag = false;
    }

    @Override
    protected void validateReporterParms() throws IOException, ParseFailureException {
        if (! this.genomeDir1.exists())
            throw new FileNotFoundException("IN-group genome source " + this.genomeDir1 + " does not exist.");
        if (! this.genomeDir2.exists())
            throw new FileNotFoundException("OUT-group genome source " + this.genomeDir2 + " does not exist.");
        if (this.min1 > 1.0)
            throw new ParseFailureException("IN-group minimum fraction cannot be greature that 1.0.");
        if (this.max2 < 0.0)
            throw new ParseFailureException("OUT-group maximum fraction cannot be less than 0.0.");
        if (this.min1 < this.max2)
            throw new ParseFailureException("Minimum IN-group fraction cannot be less than maximum OUT-group fraction.");
        // Create the classifier.
        this.classifier = this.classType.create(this);
        // Connect to the IN-group source.
        this.genomes1 = this.sourceType1.create(this.genomeDir1);
        log.info("{} genomes found in IN group.", this.genomes1.size());
        // Connect to the OUT-group source.
        this.genomes2 = this.sourceType2.create(this.genomeDir2);
        log.info("{} genomes found in OUT group.", this.genomes2.size());
    }

    @Override
    protected void runReporter(PrintWriter writer) throws Exception {
        // Create the reporting object.
        this.reporter = this.outFormat.create(writer, this);
        // Compute the class counts for the IN group.
        log.info("Computing class counts for IN group.");
        CountMap<String> counts1 = this.computeCounts(this.genomes1);
        // Compute the class counts for the OUT group.
        log.info("Computing class counts for OUT group.");
        CountMap<String> counts2 = this.computeCounts(this.genomes2);
        log.info("{} classes found for IN group, {} for OUT group.", counts1.size(), counts2.size());
        // Compute the fraction limits.
        this.reportSignatures(counts1, counts2, this.genomes1, this.genomes2);
        // Check for a BOTH request.
        if (this.bothFlag) {
            // Leave space between the two reports.
            this.reporter.space();
            this.reportSignatures(counts2, counts1, this.genomes2, this.genomes1);
        }
    }

    /**
     * Report the signatures for a particular pair of genome sets.
     *
     * @param inCounts		count map for IN group
     * @param outCounts		count map for OUT group
     * @param inSource		genome source for IN group
     * @param outSource		genome source for OUT group
     */
    private void reportSignatures(CountMap<String> inCounts, CountMap<String> outCounts, GenomeSource inSource,
            GenomeSource outSource) {
        int minCount1 = (int) Math.ceil(inSource.size() * this.min1);
        int maxCount2 = (int) Math.floor(outSource.size() * this.max2);
        // Extract all the classes that should be included in the report.
        List<String> signatures = inCounts.sortedCounts().stream().filter(x -> x.getCount() >= minCount1)
                .map(x -> x.getKey()).filter(x -> outCounts.getCount(x) <= maxCount2).collect(Collectors.toList());
        log.info("{} signature classes found.", signatures.size());
        // Create the signature name map.
        Map<String, String> sigMap = this.classifier.getNames(signatures);
        // Initialize the report.
        this.reporter.openReport(sigMap, this.genomes1.size(), this.genomes2.size());
        // Loop through the classes, creating the report.
        for (String signature : signatures)
            this.reporter.showClass(signature, inCounts.getCount(signature), outCounts.getCount(signature));
        // Finish the report.
        this.reporter.closeReport();
    }

    /**
     * @return the class counts for the specified genome group
     *
     * @param genomes	genome source for genomes to count
     */
    private CountMap<String> computeCounts(GenomeSource genomes) {
        CountMap<String> retVal = new CountMap<String>();
        for (Genome genome : genomes) {
            // Get all the classes in this genome.
            Set<String> gClasses = this.classifier.getClasses(genome);
            // Increment their counts.
            for (String gClass : gClasses)
                retVal.count(gClass);
        }
        // Return the total counts.
        return retVal;
    }

    @Override
    public File getRoleFile() {
        return this.roleFile;
    }

}
