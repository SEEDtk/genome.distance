/**
 *
 */
package org.theseed.genome.distance;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.proteins.Role;
import org.theseed.proteins.RoleMap;
import org.theseed.sequence.ProteinKmers;

/**
 * This class creates a protein kmer database for a genome's well-annotated roles.
 * It takes as input the genome itself and a RoleMap of the roles of interest.  The
 * object itself maps each role to a list of protein kmer objects for all the role's
 * proteins in the genome.
 *
 * @author Bruce Parrello
 *
 */
public class ProtMeasurer extends Measurer {

    // FIELDS
    /** map of role IDs to protein kmer lists */
    private HashMap<String, Collection<ProteinKmers>> roleKmers;
    /** number of distinct roles in this genome */
    private int roleCount;
    /** map of useful roles */
    private static RoleMap roleMap;

    /**
     * Specify the role map to use.
     *
     * @param myRoles	role map to store
     */
    public static void setRoleMap(RoleMap myRoles) {
        roleMap = myRoles;
    }

    /**
     * Construct a measurement object from a genome.
     *
     * @param genome	the genome for which a measurement object is desired
     */
    public ProtMeasurer(Genome genome) {
        super(genome);
        this.roleCount = 0;
        this.roleKmers = new HashMap<String, Collection<ProteinKmers>>();
        // Loop through the genome features.
        for (Feature feat : genome.getPegs()) {
            // Find any useful roles for this feature.
            Collection<Role> featRoles = feat.getUsefulRoles(roleMap);
            if (featRoles.size() > 0) {
                // Here the feature has at least one useful role.  Get its protein kmers.
                ProteinKmers kmerObject = new ProteinKmers(feat.getProteinTranslation());
                // Loop through the feature's roles, adding them to the hash.
                for (Role role : feat.getUsefulRoles(roleMap)) {
                    String roleId = role.getId();
                    if (! roleKmers.containsKey(roleId)) {
                        // Here we have a new role.
                        Collection<ProteinKmers> list = new ArrayList<ProteinKmers>(5);
                        list.add(kmerObject);
                        roleKmers.put(roleId, list);
                        this.roleCount++;
                    } else {
                        // Here we have an existing role, so add these kmers to the role's
                        // collection.
                        roleKmers.get(roleId).add(kmerObject);
                    }
                }
            }
        }
    }

    /**
     * Compute the percent similarity between this object's genome and another genome.
     *
     * @param genome	genome to compare to this one
     *
     * @return the percent similarity between the genomes, based on kmers for roles
     */
    @Override
    public double computePercentSimilarity(Genome genome) {
        ProtMeasurer other = new ProtMeasurer(genome);
        double retVal = computePercentSimilarity(other);
        return retVal;
    }

    /**
     * Compute the percent similarity between this object's genome and another object's genome.
     *
     * @param other		measurer for the other genome
     *
     * @return the percent similarity between the genomes, based on kmers for roles
     */
    public double computePercentSimilarity(ProtMeasurer other) {
        // This will total the best similarity for each role.
        double retVal = 0.0;
        // This will count the roles in common.
        int commonRoles = 0;
        for (String roleId : this.roleKmers.keySet()) {
            Collection<ProteinKmers> myKmerList = this.roleKmers.get(roleId);
            Collection<ProteinKmers> otherKmerList = other.roleKmers.get(roleId);
            if (otherKmerList != null) {
                // Here the role is in common and we can check the similarity.
                commonRoles++;
                double bestDistance = bestDistance(myKmerList, otherKmerList);
                retVal += 1.0 - bestDistance;
            }
        }
        // retVal is now a sum of the similarities.  Divide by the role counts and
        // convert to a percentage.
        retVal = retVal * 100.0 / (this.roleCount + other.roleCount - commonRoles);
        return retVal;
    }

    /**
     * @return the distance between two genomes based on a seed protein.
     *
     * @param other		measurer of the other genome to compare to this one
     * @param seedId	ID of the seed protein
     */
    public double computeDistance(ProtMeasurer other, String seedId) {
        Collection<ProteinKmers> myKmerList = this.roleKmers.get(seedId);
        Collection<ProteinKmers> otherKmerList = other.roleKmers.get(seedId);
        double retVal = 1.0;
        if (myKmerList != null && otherKmerList != null)
            retVal = bestDistance(myKmerList, otherKmerList);
        return retVal;
    }

    /**
     * @return the best distance between kmers in two collections
     *
     * @param myKmerList		source collection
     * @param otherKmerList		target collection
     */
    protected double bestDistance(Collection<ProteinKmers> myKmerList, Collection<ProteinKmers> otherKmerList) {
        double retVal = 1.0;
        for (ProteinKmers myKmer : myKmerList) {
            for (ProteinKmers otherKmer : otherKmerList) {
                double distance = myKmer.distance(otherKmer);
                retVal = Math.min(distance, retVal);
            }
        }
        return retVal;
    }

}