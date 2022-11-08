/**
 *
 */
package org.theseed.genome.distance.methods;

import org.theseed.genome.Genome;

/**
 * This is the base class for measurement objects that can be used to apply a single
 * measurement method to a single genome.  It has no real methods or fields, only a
 * constructor.
 *
 * @author Bruce Parrello
 *
 */
public abstract class Measurer {

    /**
     * Construct a measurement structure for a single genome.
     *
     * @param genome	genome to use
     */
    public Measurer(Genome genome) { }

}
