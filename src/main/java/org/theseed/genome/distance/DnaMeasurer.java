/**
 *
 */
package org.theseed.genome.distance;

import java.io.UnsupportedEncodingException;
import java.security.NoSuchAlgorithmException;

import org.theseed.genome.Genome;
import org.theseed.sequence.GenomeKmers;
import org.theseed.sequence.ProteinKmers;

/**
 * This is an alternate measuring method that uses contig DNA kmers to determine genome distance.
 *
 * @author Bruce Parrello
 *
 */
public class DnaMeasurer extends Measurer {

    // FIELDS
    /** genome kmer set */
    private GenomeKmers kmers;

    /**
     * Construct the DNA kmer database for a genome.
     *
     * @param genome	genome to process
     */
    public DnaMeasurer(Genome genome) {
        super(genome);
        try {
            this.kmers = new GenomeKmers(genome);
        } catch (NoSuchAlgorithmException | UnsupportedEncodingException e) {
            throw new RuntimeException(e);
        }
    }

    @Override
    public double computePercentSimilarity(Genome other) {
        GenomeKmers otherKmers;
        try {
            otherKmers = new GenomeKmers(other);
        } catch (NoSuchAlgorithmException | UnsupportedEncodingException e) {
            throw new RuntimeException(e);
        }
        int sim = this.kmers.similarity(otherKmers);
        // Convert the numerical similarity count to a percent similarity.
        double retVal = 1.0;
        if (sim == ProteinKmers.INFINITY)
            retVal = 0.0;
        else if (sim > 0.0)
            retVal = sim * 100.0 / (this.kmers.size() + otherKmers.size() - sim);
        return retVal;
    }

}
