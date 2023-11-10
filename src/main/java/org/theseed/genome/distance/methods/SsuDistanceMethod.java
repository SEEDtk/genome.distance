/**
 *
 */
package org.theseed.genome.distance.methods;

import java.util.Map;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.ParseFailureException;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.p3api.P3Genome;
import org.theseed.p3api.P3Genome.Details;
import org.theseed.sequence.RnaKmers;

/**
 * This method computes the distance between genomes using SSU kmers.  The keywords in the parameter string
 * are as follows:
 *
 * K	DNA kmer size (default 20)
 *
 * Note that even though the SSU rRNA is a DNA molecule, we use ProteinKmer because the DNA is not reversible.
 * There are usually multiple SSUs in a genome; all of them are combined into a single kmer set.
 *
 * @author Bruce Parrello
 *
 */
public class SsuDistanceMethod extends DistanceMethod {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(SsuDistanceMethod.class);
    /** DNA kmer size */
    private int kSize;

    /**
     * This class contains the data structures for a single genome.
     */
    protected class Analysis extends Measurer {

        /** sequence kmers for the SSUs */
        private RnaKmers kmers;

        /**
         * @param genome
         */
        public Analysis(Genome genome) {
            super(genome);
            // Create an blank RNA kmers object.
            this.kmers = new RnaKmers(SsuDistanceMethod.this.kSize);
            // Loop through the genome, finding SSUs.
            for (Feature feat : genome.getFeatures()) {
                if (feat.getType().contentEquals("rna") && Genome.isSSURole(feat)) {
                    // Here we have an SSU rRNA.  Fold in its sequence.
                    this.kmers.addSequence(genome.getDna(feat.getLocation()));
                }
            }
        }

    }

    @Override
    protected Measurer setupGenome(Genome genome) {
        return this.new Analysis(genome);
    }

    @Override
    protected void parseParms(Map<String, String> keywords) throws ParseFailureException {
        this.kSize = this.getIntValue(keywords, "K", 20);
    }

    @Override
    public double getDistance(Measurer measurer, Measurer other) {
        SsuDistanceMethod.Analysis m1 = (SsuDistanceMethod.Analysis) measurer;
        SsuDistanceMethod.Analysis m2 = (SsuDistanceMethod.Analysis) other;
        double retVal = m1.kmers.distance(m2.kmers);
        return retVal;
    }

    @Override
    public Details getDetailLevel() {
        return P3Genome.Details.FULL;
    }

    @Override
    public String getName() {
        return String.format("SSU_k%d", this.kSize);
    }

    @Override
    public void close() throws Exception {
    }

}
