/**
 *
 */
package org.theseed.reports;

import java.io.OutputStream;

import org.theseed.sequence.blast.BlastHit;
import org.theseed.sequence.blast.BlastHit.SeqData;

/**
 * @author Bruce Parrello
 *
 */
public class BlastHtmlReporter extends BlastReporter {

    /**
     * This enum indicates how the color of a hit is determined.
     */
    public enum ColorType {
        /** percent identity */
        ident("percent identity"),
        /** percent similarity */
        sim("percent similarity"),
        /** percent coverage of target sequence */
        covg("percent coverage");

        // FIELDS
        private String name;

        private ColorType(String name) {
            this.name = name;
        }

        /**
         * @return the color for a blast hit by the specified target sequence.
         *
         * @param target	target sequence forming the denominator
         * @param hit		blast hit whose color is needed
         */
        public Color computeColor(BlastHit.SeqData target, BlastHit hit) {
            double fraction = 1.0;
            switch (this) {
            case ident:
                fraction = ((double) hit.getNumIdentical()) / hit.getAlignLen();
                break;
            case sim:
                fraction = ((double) hit.getNumSimilar()) / hit.getAlignLen();
                break;
            case covg:
                fraction = ((double) hit.getNumSimilar()) / target.getLen();
                break;
            }
            Color retVal;
            if (fraction >= 1.0)
                retVal = Color.GRAY;
            else if (fraction >= 0.9)
                retVal = Color.DARK_GREEN.brighten((1.0 - fraction)*5);
            else if (fraction >= 0.7)
                retVal = Color.DARK_YELLOW.brighten((0.9 - fraction)*2.5);
            else if (fraction >= 0.5)
                retVal = Color.ORANGE.brighten((0.7 - fraction)*2.5);
            else
                retVal = Color.RED.brighten(0.5 - fraction);
            return retVal;
        }

        @Override
        public String toString() {
            return this.name;
        }

    }

    // FIELDS
    /** type of color display */
    private ColorType colorType;

    /**
     * Construct a report stream for HTML reports.
     *
     * @param output	output stream to contain report
     * @param sort		sequence type (query, subject) to sort on
     */
    public BlastHtmlReporter(OutputStream output, SortType sort) {
        super(output, sort);
        this.colorType = ColorType.sim;
    }

    @Override
    protected void openReport(String title) {
        // TODO start HTML report

    }

    @Override
    protected void openSection(SeqData data) {
        // TODO initialize current sort sequence

    }

    @Override
    protected void processHit(SeqData target, BlastHit hit) {
        // TODO process a BLAST hit

    }

    @Override
    protected void closeSection() {
        // TODO create the display for the current sort sequence

    }

    @Override
    protected void closeReport() {
        // TODO output the report
    }

    /**
     * Specify the scheme for computing the hit color.
     *
     * @param colorType the color computation scheme to use
     */
    public void setColorType(ColorType colorType) {
        this.colorType = colorType;
    }

}
