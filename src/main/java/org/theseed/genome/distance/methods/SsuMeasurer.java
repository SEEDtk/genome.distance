/**
 *
 */
package org.theseed.genome.distance.methods;

import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.sequence.DnaKmers;

/**
 * This measurer builds a DnaKmer hash from all of the SSUs in a genome and computes the percent of
 * kmers in common.  This works because even if the number of SSUs varies between genomes, we expect
 * the individual copies to be very similar.
 *
 * @author Bruce Parrello
 *
 */
public class SsuMeasurer extends Measurer {

    // FIELDS
    private DnaKmers ssuKmers;

    public SsuMeasurer(Genome genome) {
        super(genome);
        // Initialize the kmer hash.
        this.ssuKmers = new DnaKmers();
        // Build a DNA kmer hash from the SSU features.
        for (Feature feat : genome.getFeatures()) {
            if (feat.getType().contentEquals("rna") && Genome.isSSURole(feat)) {
                // Here we have an SSU rRNA.  Fold in its sequence.
                this.ssuKmers.addSequence(genome.getDna(feat.getLocation()));
            }
        }
    }

    @Override
    public Type getType() {
        return Type.SSU;
    }

    @Override
    public double computePercentSimilarity(Measurer otherMeasurer) {
        DnaKmers otherKmers = ((SsuMeasurer) otherMeasurer).ssuKmers;
        // Convert distance to percent similarity.
        return 100.0 - (this.ssuKmers.distance(otherKmers)) * 100.0;
    }

}
