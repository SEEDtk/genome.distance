/**
 *
 */
package org.theseed.genome.distance.methods;

import java.io.IOException;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Map;
import java.util.OptionalInt;
import java.util.stream.IntStream;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.ParseFailureException;
import org.theseed.genome.Genome;
import org.theseed.p3api.P3Genome;

/**
 * This method computes a very crude distance from the taxonomic lineage.  The 7 orthodox taxonomic levels are
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
    private static final String[] LEVELS = new String[] { "superkingdom", "phylum", "class", "order", "family", "genus", "species" };
    /** value to indicate a missing taxon ID */
    private static final int MISSING = -1;

    /**
     * This class contains the lineage data for a distance measure.
     */
    public class Analysis extends Measurer {

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

    /**
     * This is a comparator for sorting taxonomic rankings from highest to lowest.  Non-standard ranks are sorted
     * to the end in alphabetical order.
     */
    public static class RankSort implements Comparator<String> {

        @Override
        public int compare(String o1, String o2) {
            int retVal;
            int l1 = rankLevel(o1);
            int l2 = rankLevel(o2);
            if (l1 < 0) {
                if (l2 < 0)
                    retVal = o1.compareTo(o2);
                else
                    retVal = 1;
            } else if (l2 < 0)
                retVal = -1;
            else
                retVal = l2 - l1;
            return retVal;
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
        int retVal = getDiffLevel(m1, m2);
        return (retVal < 0 ? 0.0 : 1.0 / (2 << retVal));
    }

    /**
     * @return the level index for a rank, or -1 if the rank is not one of the standard ones
     *
     * @param rank	taxonomic rank name
     */
    public static int rankLevel(String rank) {
        OptionalInt idx = IntStream.range(0, LEVELS.length).filter(i -> LEVELS[i].contentEquals(rank)).findFirst();
        return (idx.orElse(MISSING));

    }

    /**
     * Compute the level at which the taxonomies of two genomes disagree.
     *
     * @param m1	analysis of first genome
     * @param m2	analysis of second genome
     *
     * @return the index of the difference level, or -1 if the taxonomies are the same species
     */
    protected int getDiffLevel(TaxonDistanceMethod.Analysis m1, TaxonDistanceMethod.Analysis m2) {
        int retVal = -1;
        for (int i = 0; retVal < 0 && i < LEVELS.length; i++) {
            if (m1.lineage[i] != m2.lineage[i])
                retVal = i;
        }
        return retVal;
    }

    /**
     * This method determines the smallest taxonomic grouping level to which two genomes belong.
     *
     * @param m1	taxonomic analysis of first genome
     * @param m2	taxonomic analysis of second genome
     *
     * @return the taxonomic point of agreement between genomes
     */
    public String getGroupingLevel(TaxonDistanceMethod.Analysis m1, TaxonDistanceMethod.Analysis m2) {
        int diffLevel = this.getDiffLevel(m1, m2);
        String retVal;
        if (diffLevel == -1) {
            // Here the genomes were the same across the whole lineage.
            retVal = LEVELS[LEVELS.length - 1];
        } else if (diffLevel == 0) {
            // Here the genomes belong to different superkingdoms.
            retVal = "cellular organism";
        } else {
            // "diffLevel" is the first point where the lineages disagree, so we back up one to get
            // the group they both belong to.
            retVal = LEVELS[diffLevel - 1];
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
