/**
 *
 */
package org.theseed.genome.distance.methods;

import java.util.List;

import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.proteins.Role;
import org.theseed.proteins.RoleMap;
import org.theseed.sequence.ProteinKmers;
import org.theseed.utils.ParseFailureException;

/**
 * This measures genome distance using a single seed protein.  Its purpose is to evaluate our standard kmer-based
 * approach to choosing representative genomes.
 *
 * @author Bruce Parrello
 *
 */
public class SeedMeasurer extends Measurer {

    // FIELDS
    /** special role map for the seed protein */
    private static RoleMap roleMap;
    /** kmer map for seed role */
    private ProteinKmers kmers;

    /**
     * Create a role map for the specified seed protein.
     *
     * @param fullMap		full role definition map
     * @param role		ID of the target seed protein
     *
     * @throws ParseFailureException
     */
    public static void setRole(RoleMap fullMap, String role) throws ParseFailureException {
        // Copy the target role to a new map.
        roleMap = new RoleMap();
        List<Role> roles = fullMap.getAllById(role);
        if (roles.isEmpty())
            throw new ParseFailureException("Role ID " + role + " not found in " + roleMap + ".");
        roleMap.putAll(roles);
    }

    /**
     * Create a measurer for the specified genome.
     *
     * @param genome
     */
    public SeedMeasurer(Genome genome) {
        super(genome);
        // Find the longest occurrence of the seed role.
        String bestProt = "";
        for (Feature peg : genome.getPegs()) {
            String prot = peg.getProteinTranslation();
            if (prot.length() > bestProt.length()) {
                // Here the protein is long enough.  Now do the slower check for a role match.
                if (! peg.getUsefulRoles(roleMap).isEmpty())
                    bestProt = prot;
            }
        }
        this.kmers = new ProteinKmers(bestProt);
    }

    @Override
    public Type getType() {
        return Measurer.Type.SEED;
    }

    @Override
    public double computePercentSimilarity(Measurer otherMeasurer) {
        ProteinKmers otherKmers = ((SeedMeasurer) otherMeasurer).kmers;
        // Get the sim score.
        double score = 1.0 - this.kmers.distance(otherKmers);
        // Convert to a percent.
        return 100.0 * score;
    }

}
