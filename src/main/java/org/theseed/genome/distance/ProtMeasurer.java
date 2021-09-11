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
     * This is a utility class to hold a count and an accumulator for creating
     * the final distance.
     */
    private static class Accumulator {

        private int count;
        private double sum;

        /**
         * Create a blank accumulator.
         */
        protected Accumulator() {
            this.count = 0;
            this.sum = 0.0;
        }

        /**
         * Add a value to the sum
         *
         * @param value		value to add
         */
        protected void add(double value) {
            this.sum += value;
            this.count++;
        }

        /**
         * Merge two accumulators.
         *
         * @param other		other accumulator to merge
         */
        protected void merge(Accumulator other) {
            this.sum += other.sum;
            this.count += other.count;
        }

        /**
         * Convert this from a sum of distances to a percent similarity.
         *
         * @param count1	number of roles in the first genome
         * @param count2	number of roles in the second genome
         */
        protected double percent(int count1, int count2) {
            return (this.count - this.sum) * 100.0 / (count1 + count2 - this.count);
        }
    }

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
        // Process all the roles.
        Accumulator accum = this.roleKmers.keySet().parallelStream()
                .collect(Accumulator::new, (x, r) -> this.processRole(other, x, r),
                        (x1,x2) -> x1.merge(x2));
        double retVal = accum.percent(this.roleCount, other.roleCount);
        return retVal;
    }

    /**
     * Process a single role for the comparison.
     *
     * @param other		measurer for the other genome
     * @param accum		accumulator in which to accumulate the results
     * @param roleId	ID of the role to process
     */
    private void processRole(ProtMeasurer other, Accumulator accum, String roleId) {
        Collection<ProteinKmers> myKmerList = this.roleKmers.get(roleId);
        Collection<ProteinKmers> otherKmerList = other.roleKmers.get(roleId);
        if (otherKmerList != null) {
            // Here the role is in common and we can check the similarity.
            double bestDistance = this.bestDistance(myKmerList, otherKmerList);
            accum.add(bestDistance);
        }
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
            retVal = this.bestDistance(myKmerList, otherKmerList);
        return retVal;
    }

    /**
     * Compute the lowest distance between items in two kmer lists.  This is not parallelized
     * because we expect the two lists to be very small, usually one or two kmer sets.
     *
     * @param myKmerList		source collection
     * @param otherKmerList		target collection
     *
     * @return the best distance between kmers in two collections
     */
    protected double bestDistance(Collection<ProteinKmers> myKmerList,
            Collection<ProteinKmers> otherKmerList) {
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
