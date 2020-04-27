/**
 *
 */
package org.theseed.reports;

import java.io.Closeable;
import java.io.OutputStream;
import java.io.PrintWriter;
import org.theseed.genome.Genome;
import org.theseed.sequence.blast.BlastHit;

/**
 * This is the base class for reports about feature-to-contig matches.
 *
 * @author Bruce Parrello
 *
 */
public abstract class MatchReporter implements Closeable, AutoCloseable {

    // FIELDS
    /** genome containing the features */
    private Genome fidGenome;
    /** ID of the current contig */
    private String contigId;
    /** length of the current contig */
    private int contigLen;
    /** output stream for the report */
    private PrintWriter writer;

    /**
     * Enumeration for type of report.
     */
    public enum Type {
        TABLE,
        HTML;
    }

    /**
     * Create a new match reporter.
     *
     * @param type		type of report
     * @param genome	genome containing the matching features
     * @param output	output stream to contain report
     */
    public static MatchReporter create(Type type, Genome genome, OutputStream output) {
        MatchReporter retVal = null;
        switch (type) {
        case TABLE:
            retVal = new MatchTableReporter(genome, output);
            break;
        case HTML:
            retVal = new MatchHtmlReporter(genome, output);
            break;
        }
        return retVal;
    }

    /**
     * Construct a blank report.
     *
     * @param genome	genome containing the matching features
     */
    protected MatchReporter(Genome genome, OutputStream output) {
        this.fidGenome = genome;
        this.writer = new PrintWriter(output);
        // Denote we don't have a contig yet.
        this.contigId = null;
        // Start the report.
        this.openReport();
    }

    /**
     * Begin the report.
     */
    protected abstract void openReport();

    /**
     * Specify the current contig for this section of the report.
     *
     * @param contigId	ID of the contig
     * @param length	length of the contig
     */
    public void startContig(String contigId, int length) {
        // Clean up the old contig (if any).
        if (this.contigId != null)
            this.closeContig();
        // Store the new contig's information.
        this.contigId = contigId;
        this.contigLen = length;
        // Start the new contig.
        this.openContig();
    }

    /**
     * Begin reporting on a new contig.
     */
    protected abstract void openContig();

    /**
     * Specify a BLAST hit for the current contig.
     *
     * @param hit	BLAST hit specification
     */
    public void recordHit(BlastHit hit) {
        this.processHit(hit);
    }

    /**
     * Process a BLAST hit.
     *
     * @param hit	descriptor representing the BLAST hit
     */
    protected abstract void processHit(BlastHit hit);

    /**
     * Finish reporting of the current contig.
     */
    protected abstract void closeContig();

    /**
     * Finish the entire report.
     */
    protected abstract void closeReport();

    @Override
    public void close() {
        // Make sure we've terminated the current contig.
        if (contigId != null)
            this.closeContig();
        // Terminate the whole report.
        this.closeReport();
        // Close the output stream.
        this.writer.close();
    }

    /**
     * @return the genome containing the features on this report
     */
    public Genome getGenome() {
        return fidGenome;
    }

    /**
     * @return the ID of the current contig
     */
    public String getContigId() {
        return contigId;
    }

    /**
     * @return the length of the current contig
     */
    public int getContigLen() {
        return contigLen;
    }

    /**
     * Write a formatted output line.
     */
    protected void print(String format, Object... args) {
        this.writer.format(format, args);
    }

    /**
     * Write an unformatted output line.
     */
    protected void println(String line) {
        this.writer.println(line);
    }


}
