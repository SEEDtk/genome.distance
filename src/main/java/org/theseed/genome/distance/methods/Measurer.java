/**
 *
 */
package org.theseed.genome.distance.methods;

import org.theseed.genome.Genome;

/**
 * This is the base class for measurement objects that can be used to apply a single
 * measurement method to a single genome.
 *
 * @author Bruce Parrello
 *
 */
public abstract class Measurer {

    // FIELDS
    /** ID of the genome being measured */
    private String genomeId;

    /**
     * Construct a measurement structure for a single genome.
     *
     * @param genome	genome to use
     */
    public Measurer(Genome genome) {
        this.genomeId = genome.getId();
    }

    /**
     * @return the ID of the genome being measured
     */
    public String getGenomeId() {
        return this.genomeId;
    }

}
