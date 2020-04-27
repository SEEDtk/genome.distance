/**
 *
 */
package org.theseed.reports;

import java.io.OutputStream;

import org.theseed.genome.Genome;
import org.theseed.locations.Location;
import org.theseed.sequence.blast.BlastHit;

/**
 * @author parrello
 *
 */
public class MatchTableReporter extends MatchReporter {

    protected MatchTableReporter(Genome genome, OutputStream output) {
        super(genome, output);
    }

    @Override
    protected void openReport() {
        // Print the output header.
        this.println("seq_id\tseq_len\tstrand\tleft\tright\tpeg_id\tpct_coverage\te_value\talign_len\tfunction");

    }

    @Override
    protected void openContig() {
        // Space before this contig's hits.
        System.out.println();
    }

    @Override
    protected void processHit(BlastHit hit) {
        // Get the location in the DNA sequence of the hit.
        Location qLoc = hit.getQueryLoc();
        // Write the output line.
        System.out.format("%s\t%d\t%c\t%d\t%d\t%s\t%4.2f\t%4.3e\t%d\t%s%n",
                hit.getQueryId(), hit.getQueryLen(), qLoc.getDir(), qLoc.getLeft(), qLoc.getRight(),
                hit.getSubjectId(), hit.getSubjectPercentMatch(), hit.getEvalue(), hit.getAlignLen(),
                hit.getSubjectDef());
    }

    @Override
    protected void closeContig() {
    }

    @Override
    protected void closeReport() {
    }

}
