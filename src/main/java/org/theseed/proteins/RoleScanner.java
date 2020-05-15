/**
 *
 */
package org.theseed.proteins;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.genome.GenomeDirectory;

/**
 * This package maintains a role map that can be created from a genome or directory of genomes.
 * Each genome is scanned and all non-hypothetical roles are added to the map.
 *
 * @author Bruce Parrello
 *
 */
public class RoleScanner extends RoleMap {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(RoleScanner.class);

    /**
     * Add a genome's roles to this map.
     *
     * @param gto	genome whose roles are to be added to this map
     */
    public void addGenome(Genome gto) {
        // Loop through the proteins in this genome.
        for (Feature peg : gto.getPegs()) {
            // Get all of the protein's roles and insure they are in the genome.
            for (String roleDesc : peg.getRoles())
                this.findOrInsert(roleDesc);
        }
        log.info("{} roles in map after scan of {}.", this.fullSize(), gto);
    }

    /**
     * Add the roles for all the genomes in the specified directory to this map.
     *
     * @param gDir	directory of genome to scan
     */
    public void addGenomes(GenomeDirectory gDir) {
        for (Genome gto : gDir)
            this.addGenome(gto);
    }

}
