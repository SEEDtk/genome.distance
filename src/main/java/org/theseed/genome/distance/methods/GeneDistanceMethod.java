/**
 *
 */
package org.theseed.genome.distance.methods;

import org.theseed.genome.Feature;
import org.theseed.p3api.P3Genome;

/**
 * This is a sequence distance method that uses the DNA sequence.  It is an attempt to compute true ANI distance.
 *
 * The keyword parameters are as follows:
 *
 * K		kmer size (default is 20)
 * roles	comma-delimited list of role IDs (default is to use all)
 * penalty	"Y" to penalize for a role that is present in only one genome, "N" to ignore (default is "N")
 *
 * @author Bruce Parrello
 *
 */
public class GeneDistanceMethod extends SequenceDistanceMethod {

    @Override
    protected String getSequence(Feature feat) {
        return feat.getDna();
    }

    @Override
    protected boolean getPenaltyDefault() {
        return false;
    }

    @Override
    protected int getDefaultKmerSize() {
        return 20;
    }

    @Override
    public P3Genome.Details getDetailLevel() {
        return P3Genome.Details.FULL;
    }

    @Override
    public String getName() {
        String retVal = String.format("Genes_k%d:%s", this.kSize, this.roleSource);
        if (this.penalty) retVal += ",w/penalty";
        return retVal;
    }

}
