package org.theseed.genome.distance;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;

import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.BaseProcessor;
import org.theseed.io.TabbedLineReader;
import org.theseed.sequence.ProteinKmers;

public abstract class ProteinKmerReader extends BaseProcessor {

    // FIELDS

    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(WidthProcessor.class);
    /** input stream */
    private TabbedLineReader inStream;
    /** group ID column */
    private int idIdx;
    /** protein sequence column */
    private int protIdx;

    // COMMAND-LINE OPTIONS

    /** kmer length */
    @Option(name = "-K", aliases = { "--kmer", "--kmerSize" }, metaVar = "12", usage = "protein kmer size")
    private int kmerSize;
    /** input file */
    @Option(name = "-i", aliases = { "--input" }, metaVar = "families.tbl", usage = "input file (if not STDIN)")
    private File inFile;
    /** group ID column name */
    @Option(name = "-c", aliases = { "--col", "--groupCol" }, metaVar = "pgfam_id", usage = "group ID column index (1-based) or name")
    private String idColumn;
    /** protein sequence column name */
    @Option(name = "-p", aliases = { "--prot", "--protCol" }, metaVar = "0", usage = "protein sequence column index (1-based) or name")
    private String protColumn;

    public ProteinKmerReader() {
        super();
    }

    /**
     * Initialize the parameters for the protein input file.
     */
    protected void initProteinParms() {
        this.kmerSize = 8;
        this.inFile = null;
        this.idColumn = "1";
        this.protColumn = "aa_sequence";
    }

    /**
     * Validate the parameters for the protein input.
     *
     * @throws IOException
     * @throws FileNotFoundException
     */
    protected void validateProteinParms() throws IOException, FileNotFoundException {
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
        // Set the kmer size.
        ProteinKmers.setKmerSize(this.kmerSize);
    }

    /**
     * @return the kmer object for the protein on the current line
     *
     * @param line	current input line
     */
    protected ProteinKmers getProteinKmers(TabbedLineReader.Line line) {
        return new ProteinKmers(line.get(protIdx));
    }

    /**
     * @return the group ID for the protein on the current line
     *
     * @param line	current input line
     */
    protected String getGroupId(TabbedLineReader.Line line) {
        return line.get(this.idIdx);
    }

    /**
     * @return the input stream
     */
    protected TabbedLineReader input() {
        return this.inStream;
    }

    /**
     * Process the proteins.
     */
    protected abstract void processProteins() throws IOException;

    @Override
    public void runCommand() throws Exception {
        try {
            this.processProteins();
        } finally {
            this.inStream.close();
        }
    }
}
