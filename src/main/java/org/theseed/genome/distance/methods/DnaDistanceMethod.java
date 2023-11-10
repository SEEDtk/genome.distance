/**
 *
 */
package org.theseed.genome.distance.methods;

import java.io.UnsupportedEncodingException;
import java.security.NoSuchAlgorithmException;
import java.util.Map;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.ParseFailureException;
import org.theseed.genome.Genome;
import org.theseed.p3api.P3Genome;
import org.theseed.p3api.P3Genome.Details;
import org.theseed.sequence.GenomeKmers;

/**
 * This method computes the distance between genomes using whole-DNA kmer comparison.  The keywords in the
 * parameter string are as follows.
 *
 * K	DNA kmer size (default 20)
 *
 * @author Bruce Parrello
 *
 */
public class DnaDistanceMethod extends DistanceMethod {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(DnaDistanceMethod.class);
    /** DNA kmer size */
    private int kSize;

    /**
     * This class contains the measurement structures for a single genome.
     */
    protected class Analysis extends Measurer {

        /** genome kmer set */
        private GenomeKmers kmers;

        /**
         * Analyze the specified genome.
         *
         * @param genome	genome to parse
         */
        protected Analysis(Genome genome) {
            super(genome);
            try {
                this.kmers = new GenomeKmers(genome, DnaDistanceMethod.this.kSize);
            } catch (NoSuchAlgorithmException | UnsupportedEncodingException e) {
                throw new RuntimeException(e);
            }
        }

    }


    @Override
    protected Measurer setupGenome(Genome genome) {
        return this.new Analysis(genome);
    }

    @Override
    protected void parseParms(Map<String, String> keywords) throws ParseFailureException {
        this.kSize = super.getIntValue(keywords, "K", 20);
    }

    @Override
    public double getDistance(Measurer measurer, Measurer other) {
        DnaDistanceMethod.Analysis m1 = (DnaDistanceMethod.Analysis) measurer;
        DnaDistanceMethod.Analysis m2 = (DnaDistanceMethod.Analysis) other;
        double retVal = m1.kmers.distance(m2.kmers);
        return retVal;
    }

    @Override
    public Details getDetailLevel() {
        return P3Genome.Details.CONTIGS;
    }

    @Override
    public String getName() {
        return String.format("Contigs_k%d", this.kSize);
    }

    @Override
    public int hashCode() {
        final int prime = 31;
        int result = 1;
        result = prime * result + this.kSize;
        return result;
    }

    @Override
    public boolean equals(Object obj) {
        if (this == obj) {
            return true;
        }
        if (!(obj instanceof DnaDistanceMethod)) {
            return false;
        }
        DnaDistanceMethod other = (DnaDistanceMethod) obj;
        if (this.kSize != other.kSize) {
            return false;
        }
        return true;
    }

    @Override
    public void close() throws Exception {
    }

}
