/**
 *
 */
package org.theseed.genome.distance.methods;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.proteins.RoleMap;

/**
 * This object represents a measurement method, and is used for commands that apply multiple methods.
 *
 * @author Bruce Parrello
 *
 */
public class MethodDescriptor implements Measurer.IParms {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(MethodDescriptor.class);
    /** comparison method type */
    private Measurer.Type type;
    /** seed protein ID */
    private String seedId;
    /** role definition map */
    private static RoleMap roleMap;
    /** DNA kmer size */
    private static int dnaKmerSize;

    @Override
    public int getKmerSize() {
        // TODO code for getKmerSize
        return 0;
    }

    @Override
    public String getSeedId() {
        // TODO code for getSeedId
        return null;
    }

    @Override
    public RoleMap getRoleMap() {
        // TODO code for getRoleMap
        return null;
    }


    // TODO constructors and methods for MethodDescriptor
}
