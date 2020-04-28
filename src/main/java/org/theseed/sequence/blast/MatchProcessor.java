/**
 *
 */
package org.theseed.sequence.blast;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.apache.commons.io.FileUtils;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Genome;
import org.theseed.reports.MatchReporter;
import org.theseed.sequence.DnaInputStream;
import org.theseed.sequence.SequenceDataStream;
import org.theseed.sequence.blast.BlastDB;
import org.theseed.sequence.blast.BlastHit;
import org.theseed.sequence.blast.BlastParms;
import org.theseed.sequence.blast.DnaBlastDB;
import org.theseed.sequence.blast.ProteinBlastDB;
import org.theseed.utils.BaseProcessor;

/**
 * This command takes DNA from the input and BLASTS it against the pegs of a genome to determine which protein
 * functions are present in each incoming sequence.  The intent is to determine what parts of the genome are
 * covered by the incoming sequences.
 *
 * The positional parameter is the name of the genome file.  The DNA FASTA file should come in via the standard
 * input.
 *
 * The command-line options are as follows.
 *
 * -v	show more detail on the log
 * -h	display usage information
 * -i	name of the input file for the DNA (the default is STDIN)
 * -b	number of input sequences to process at a time (the default is 10)
 *
 * --maxE		maximum permissible e-value (the default is 1e-10)
 * --minPct		minimum percent coverage of PEG for a match (the default is 75)
 * --tempDir	temporary directory for BLAST databases; the default is "Temp" in the current directory
 * --format		format for the output report (TABLE or HTML)
 * --dna		use DNA sequences of genome features instead of protein sequences of genome PEGs
 *
 * @author Bruce Parrello
 *
 */
public class MatchProcessor extends BaseProcessor {

    // FIELDS
    /** logging facility */
    protected static final Logger log = LoggerFactory.getLogger(MatchProcessor.class);
    /** DNA input stream */
    private DnaInputStream inStream;
    /** output reporter */
    private MatchReporter reporter;
    /** input genome */
    private Genome genome;

    // COMMAND-LINE OPTION

    /** input file */
    @Option(name = "-i", aliases = { "--input" }, metaVar = "rna.fasta", usage = "input file (if not STDIN)")
    private File inFile;

    /** number of DNA sequences to submit in each BLAST call */
    @Option(name = "-b", aliases = { "--batchSize", "--batch" }, metaVar = "1",
            usage = "number of input sequences to submit to each BLAST call")
    private int batchSize;

    /** maximum permissible e-value */
    @Option(name = "--maxE", aliases = { "--evalue" }, metaVar = "1e-20", usage = "maximum permissible e-value for a match")
    private double eValue;

    /** minimum permissible PEG coverage */
    @Option(name = "--minPct", aliases = { "--covg", "--coverage" }, metaVar = "67.7",
            usage = "minimum permissible coverage of the target proteins")
    private double minCoverage;

    /** temporary directory for BLAST database */
    @Option(name = "--tempDir", metaVar = "Tmp", usage = "temporary directory for BLAST databases")
    private File tempDir;

    /** output report type */
    @Option(name = "--format", usage = "output format for the report")
    private MatchReporter.Type outputFormat;

    /** genome processing mode */
    @Option(name = "--dna", usage = "if specified, the genome feature DNA will be used, instead of the protein sequences")
    private boolean dnaMode;

    /** file containing the target genome */
    @Argument(index = 0, metaVar = "genome.gto", usage = "target genome file containing the proteins",
            required = true)
    private File genomeFile;

    @Override
    protected void setDefaults() {
        this.inFile = null;
        this.batchSize = 10;
        this.eValue = 1e-10;
        this.minCoverage = 75.0;
        this.tempDir = new File(System.getProperty("user.dir"), "Temp");
        this.outputFormat = MatchReporter.Type.HTML;
        this.dnaMode = false;
    }

    @Override
    protected boolean validateParms() throws IOException {
        // Read the genome.
        log.info("Loading genome from {}.", this.genomeFile);
        this.genome = new Genome(this.genomeFile);
        // Connect to the DNA input stream.
        if (this.inFile == null) {
            this.inStream = new DnaInputStream(System.in, this.genome.getGeneticCode());
            log.info("DNA sequences will be read from standard input.");
        } else {
            this.inStream = new DnaInputStream(this.inFile, this.genome.getGeneticCode());
            log.info("DNA sequences will be read from {}.", this.inFile);
        }
        if (! this.tempDir.isDirectory()) {
            log.info("Creating temporary file directory {}.", this.tempDir);
            FileUtils.forceMkdir(this.tempDir);
        }
        // Validate the parameters.
        if (this.eValue >= 1.0)
            throw new IllegalArgumentException("Invalid eValue specified.  Must be less than 1.");
        if (this.minCoverage > 100.0)
            throw new IllegalArgumentException("Minimum coverage percent must be 100 or less.");
        if (this.batchSize <= 0)
            throw new IllegalArgumentException("Batch size must be 1 or more.");
        if (! this.genomeFile.canRead())
            throw new FileNotFoundException("Input file" + this.genomeFile + " is not found or unreadable.");
        return true;
    }

    @Override
    public void run() {
        try {
            // Create the BLAST database.
            BlastDB blastDB = this.createBlastDb();
            // Create the BLAST parameters.
            BlastParms parms = new BlastParms().maxE(this.eValue).pctLenOfSubject(this.minCoverage);
            // Now we loop through the input FASTA stream, building batches.
            int batchCount = 0;
            int seqCount = 0;
            Iterator<SequenceDataStream> batcher = this.inStream.batchIterator(batchSize);
            while (batcher.hasNext()) {
                SequenceDataStream batch = batcher.next();
                batchCount++;
                seqCount += batch.size();
                log.info("Processing input batch {}.", batchCount);
                this.processBatch(batch, blastDB, parms);
            }
            log.info("All done. {} sequences in {} batches processed.", seqCount, batchCount);
            reporter.close();
        } catch (Exception e) {
            e.printStackTrace(System.err);
        } finally {
            this.inStream.close();
        }
    }

    /**
     * Process a batch of input sequences and write their output.
     *
     * @param batch		batch of DNA sequences to process
     * @param blastDB	blast database to blast against
     * @param parms		blast parameters
     *
     * @throws InterruptedException
     * @throws IOException
     */
    private void processBatch(SequenceDataStream batch, BlastDB blastDB, BlastParms parms) throws IOException, InterruptedException {
        // Perform the BLAST.
        Map<String, List<BlastHit>> hitMap = BlastHit.sort(batch.blast(blastDB, parms));
        // Get a sorted list of the query sequence IDs that produced results.
        ArrayList<String> queryList = new ArrayList<String>(hitMap.keySet());
        queryList.sort(null);
        // Create a comparator for sorting the blast hits by location.
        Comparator<BlastHit> sortByLoc = new BlastHit.ByLoc(BlastHit.QUERY);
        for (String seqId : queryList) {
            // Get the hits for this sequence and sort them by location.
            List<BlastHit> hits = hitMap.get(seqId);
            hits.sort(sortByLoc);
            this.reporter.startContig(seqId, hits.get(0).getQueryLen());
            // Output the hits.
            for (BlastHit hit : hits)
                this.reporter.recordHit(hit);
        }
    }

    /**
     * Read in the genome and create the blast database.
     *
     * @throws IOException
     * @throws InterruptedException
     */
    private BlastDB createBlastDb() throws IOException, InterruptedException {
        this.reporter = MatchReporter.create(this.outputFormat, genome, System.out);
        File tempFile = File.createTempFile("blast", ".fasta", this.tempDir);
        BlastDB retVal;
        if (this.dnaMode)
            retVal = DnaBlastDB.createFromFeatures(tempFile, genome);
        else
            retVal = ProteinBlastDB.create(tempFile, genome);
        retVal.deleteOnExit();
        return retVal;
    }

}
