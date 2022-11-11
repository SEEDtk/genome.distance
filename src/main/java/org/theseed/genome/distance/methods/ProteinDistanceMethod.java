/**
 *
 */
package org.theseed.genome.distance.methods;

import org.theseed.genome.Feature;
import org.theseed.p3api.P3Genome;

/**
 * This is a sequence method that uses protein distance.  It is the original whole-genome distance method and
 * is designed to approximate the functional evolution of the organism.
 *
 * The keyword parameters are as follows:
 *
 * K		kmer size (default is 8)
 * roles	comma-delimited list of role IDs (default is to use all)
 * penalty	"Y" to penalize for a role that is present in only one genome, "N" to ignore (default is "Y")
 *
 * @author Bruce Parrello
 *
 */
public class ProteinDistanceMethod extends SequenceDistanceMethod {

    @Override
    protected String getSequence(Feature feat) {
        return feat.getProteinTranslation();
    }

    @Override
    protected boolean getPenaltyDefault() {
        return true;
    }

    @Override
    protected int getDefaultKmerSize() {
        return 8;
    }

    @Override
    public P3Genome.Details getDetailLevel() {
        return P3Genome.Details.PROTEINS;
    }

    @Override
    public String getName() {
        String retVal = String.format("Proteins_k%d:%s", this.kSize, this.roleSource);
        if (! this.penalty) retVal += ",noPenalty";
        return retVal;
    }

}
