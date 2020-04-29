/**
 *
 */
package org.theseed.reports;

import java.io.OutputStream;

import org.theseed.locations.Location;
import org.theseed.sequence.blast.BlastHit;
import org.theseed.sequence.blast.BlastHit.SeqData;

/**
 * This is a tabular report that shows the blast results in alignment format.  The first column gives the
 * eValue, and the remaining columns show data from the anchor sequence atop the same data from the
 * target sequence.  Only the longest alignment for each query sequence is kept.
 *
 * @author Bruce Parrello
 *
 */
public class BlastAlignReporter extends BlastReporter {

    // FIELDS
    /** best hit for this anchor */
    private BlastHit bestHit;

    /**
     * Construct an alignment report.
     *
     * @param output	output stream to receive the report
     * @param sort		sort type (query or subject)
     */
    public BlastAlignReporter(OutputStream output, SortType sort) {
        super(output, sort);
    }

    @Override
    protected void openReport(String title) {
        // Display the headings.
        this.println("e_value\tid\tdescription\tstart\tend\tmatch_len\ttot_len\tcoverage\tsequence");
    }

    @Override
    protected void openSection(SeqData data) {
        this.bestHit = null;
    }

    @Override
    protected void processHit(SeqData target, SeqData anchor, BlastHit hit) {
        if (this.bestHit == null) {
            this.bestHit = hit;
        } else if (this.bestHit.getNumSimilar() < hit.getNumSimilar()) {
            this.bestHit = hit;
        }
    }

    /**
     * Display an output line.
     *
     * @param eColumn	text for the e-value column
     * @param data		sequence data for this line
     * @param hit		blast hit being displayed
     */
    private void showData(String eColumn, SeqData data, BlastHit hit) {
        Location loc = data.getLoc();
        double coverage = hit.getNumSimilar() * 100.0 / data.getLen();
        this.print("%s\t%s\t%s\t%d\t%d\t%d\t%d\t%4.1f\t%s", eColumn, data.getId(), data.getDef(),
                loc.getBegin(), loc.getEnd(), loc.getLength(), data.getLen(),
                coverage, data.getAlignment());
    }

    @Override
    protected void closeSection() {
        if (this.bestHit != null) {
            // Space before this group.
            this.println();
            // Write the e-value and the anchor sequence.
            this.showData(String.format("%4.2e", this.bestHit.getEvalue()),
                    this.getSortType().data(this.bestHit), this.bestHit);
            // Write the target sequence with the e-column blank.
            this.showData("", this.getSortType().target(this.bestHit), this.bestHit);
        }
    }

    @Override
    protected void closeReport() {
    }

}
