/**
 *
 */
package org.theseed.reports;

import java.io.OutputStream;
import java.util.ArrayList;
import java.util.List;

import org.theseed.genome.Genome;
import org.theseed.locations.Location;
import org.theseed.sequence.blast.BlastHit;

import j2html.tags.ContainerTag;
import j2html.tags.DomContent;
import static j2html.TagCreator.*;

/**
 * This produces a visual report of the matches against a genome's pegs.
 *
 * @author Bruce Parrello
 *
 */
public class MatchHtmlReporter extends MatchReporter {

    // FIELDS
    /** current contig object */
    private HtmlFullSequence contig;
    /** list of contig displays */
    private List<DomContent> pieces;

    /**
     * Construct a new HTML match report.
     *
     * @param genome	target genome
     * @param output	output stream
     */
    public MatchHtmlReporter(Genome genome, OutputStream output) {
        super(genome, output);
    }

    @Override
    protected void openReport() {
        // We will store the parts of the report in here.
        this.pieces = new ArrayList<DomContent>();
        // Create a title and put it on the pieces list.
        DomContent title = h1(String.format("Match Report for %s", this.getGenome().toString()));
        this.pieces.add(title);
    }

    @Override
    protected void openContig() {
        // Create the HTML contig for this new contig.
        this.contig = new HtmlFullSequence(1, this.getContigLen(), this.getContigId());
    }

    @Override
    protected void processHit(BlastHit hit) {
        String label = String.format("ident=%3.2f%% E=%4.3e matched=%d/%d %s", hit.getPercentIdentity(),
                hit.getEvalue(), hit.getNumSimilar(), hit.getSubjectLen(), hit.getSubjectDef());
        Color color = MatchHtmlReporter.computeColor(hit.getSubjectPercentMatch());
        Location qLoc = hit.getQueryLoc();
        char dir = (qLoc.getDir() != hit.getSubjectLoc().getDir() ? '-' : '+');
        HtmlHitSequence hFeat = new HtmlHitSequence(hit.getSubjectId(), label,
                qLoc, dir, color);
        this.contig.add(hFeat);
    }

    /**
     * @return the color for the specified similarity level
     *
     * @param percentSimilarity	input similarity level
     */
    private static Color computeColor(double percentSimilarity) {
        // The default is gray if the similarity is over 100 (which is an error).
        double factor = 0.0;
        Color base = Color.GRAY;
        if (percentSimilarity <= 100.0) {
            if (percentSimilarity > 90.0) {
                base = Color.DARK_GREEN;
                factor = (100.0 - percentSimilarity) / 20.0;
            } else if (percentSimilarity > 70) {
                base = Color.DARK_YELLOW;
                factor = (90.0 - percentSimilarity) / 40.0;
            } else if (percentSimilarity > 50) {
                base = Color.ORANGE;
                factor = (70.0 - percentSimilarity) / 40.0;
            } else {
                base = Color.RED;
                factor = (50.0 - percentSimilarity) / 100.0;
            }
        }
        return base.brighten(factor);
    }

    @Override
    protected void closeContig() {
        DomContent contigHtml = this.contig.draw(this.getGenome().getLinker());
        this.pieces.add(contigHtml);
    }

    @Override
    protected void closeReport() {
        ContainerTag page = html(head(title(this.getGenome().getName())),
                body().with(this.pieces));
        this.println(page.render());
    }

}
