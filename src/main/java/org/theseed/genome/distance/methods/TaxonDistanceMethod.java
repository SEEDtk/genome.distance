/**
 *
 */
package org.theseed.genome.distance.methods;

import java.io.IOException;
import java.util.Arrays;
import java.util.Map;
import java.util.OptionalInt;
import java.util.stream.IntStream;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Genome;
import org.theseed.p3api.P3Genome;
import org.theseed.utils.ParseFailureException;

/**
 * This method computes a very crude distance from the taxonomic lineage.  The 9 orthodox taxonomic levels are
 * examined one by one.  The process stops on the first mismatch.  A mismatch on superkingdom is a distance
 * of 1.0.  A mismatch on a lower level is worth half the value of the previous level.
 *
 * There are no keywords.
 *
 * @author Bruce Parrello
 *
 */
public class TaxonDistanceMethod extends DistanceMethod {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(TaxonDistanceMethod.class);
    /** ranking list */
    private static final String[] LEVELS = new String[] { "superkingdom", "phylum", "class", "order", "family", "genus", "species", "strain" };
    /** value to indicate a missing taxon ID */
    private static final int MISSING = -1;

    /**
     * This class contains the lineage data for a distance measure.
     */
    protected class Analysis extends Measurer {

        /** array of taxon IDs at each level */
        private int[] lineage;

        /**
         * Create a taxonomic comparator for the specified genome.
         *
         * @param genome	genome to parse
         */
        public Analysis(Genome genome) {
            super(genome);
            // Create an empty lineage array.
            this.lineage = new int[LEVELS.length];
            Arrays.fill(this.lineage, MISSING);
            // Iterate through the lineage, storing tax IDs.
            var iter = genome.taxonomy();
            while (iter.hasNext()) {
                var taxItem = iter.next();
                final String rank = taxItem.getRank();
                OptionalInt idx = IntStream.range(0, LEVELS.length).filter(i -> LEVELS[i].contentEquals(rank)).findFirst();
                idx.ifPresent(i -> lineage[i] = taxItem.getId());
            }
        }

    }

    @Override
    protected Measurer setupGenome(Genome genome) {
        return this.new Analysis(genome);
    }

    @Override
    protected void parseParms(Map<String, String> keywords) throws ParseFailureException, IOException {
    }

    @Override
    public double getDistance(Measurer measurer, Measurer other) {
        TaxonDistanceMethod.Analysis m1 = (TaxonDistanceMethod.Analysis) measurer;
        TaxonDistanceMethod.Analysis m2 = (TaxonDistanceMethod.Analysis) other;
        double retVal = 0.0;
        for (int i = 0; retVal == 0 && i < LEVELS.length; i++) {
            if (m1.lineage[i] != m2.lineage[i])
                retVal = 1.0 / (2 << i);
        }
        return retVal;
    }

    @Override
    public P3Genome.Details getDetailLevel() {
        return P3Genome.Details.STRUCTURE_ONLY;
    }

    @Override
    public String getName() {
        return "Taxonomic";
    }

    @Override
    public void close() throws Exception {
    }

}
