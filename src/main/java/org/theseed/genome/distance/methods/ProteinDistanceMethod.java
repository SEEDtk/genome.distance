/**
 *
 */
package org.theseed.genome.distance.methods;

import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

import org.apache.commons.lang3.StringUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.p3api.P3Genome;
import org.theseed.p3api.P3Genome.Details;
import org.theseed.proteins.Role;
import org.theseed.proteins.RoleMap;
import org.theseed.sequence.ProteinKmers;
import org.theseed.utils.ParseFailureException;

/**
 * This distance method does a kmer comparison of identical roles in genomes.  For each role, the longest protein
 * is kept and converted into protein kmers.  This gives us greater precision than ProfileDistanceMethod but removes
 * the influence of non-coding regions.  Note that the use of the longest protein simplifies the algorithm, but it
 * leaves us vulnerable to contamination:  a long protein from an alien species will appear very different from what
 * is supposed to be there.
 *
 * The default for this method is to use all the roles in the standard map; however, a small list of roles can be
 * specified instead.  This allows a single seed protein role (as is used in repgen sets) or a group of common
 * universal roles (as is used in hammer generation).
 *
 * The keyword parameters are as follows:
 *
 * K		protein kmer size (default 8)
 * roles	comma-delimited list of role IDs (default is to use all)
 *
 * @author Bruce Parrello
 *
 */
public class ProteinDistanceMethod extends DistanceMethod {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(ProteinDistanceMethod.class);
    /** role map to use */
    private RoleMap protRoles;
    /** protein kmer size */
    private int kSize;
    /** name for role source */
    private String roleSource;

    /**
     * This object holds the data structures for an individual genome.
     */
    protected class Analysis extends Measurer {

        /** hash of role IDs to kmers */
        private Map<String, ProteinKmers> kmerMap;

        /**
         * Analyze a genome's proteins to produce the kmers.
         *
         * @param genome
         */
        public Analysis(Genome genome) {
            super(genome);
            var protMap = new ConcurrentHashMap<String, String>();
            // Loop through the genome features, scanning for useful pegs and storing the kmers.
            genome.getFeatures().parallelStream().forEach(x -> this.scanFeature(protMap, x));
            // Roll up the proteins into a protein kmer map.
            this.kmerMap = new HashMap<String, ProteinKmers>(protMap.size() * 4 / 3);
            for (var protEntry : protMap.entrySet()) {
                ProteinKmers kmers = new ProteinKmers(protEntry.getValue(), ProteinDistanceMethod.this.kSize);
                this.kmerMap.put(protEntry.getKey(), kmers);
            }
        }

        /**
         * Analyze a feature, and if it is a peg with useful roles, store its kmers in the map.
         *
         * @param protMap	target map of role IDs to protein strings
         * @param feat		feature to analyze
         */
        private void scanFeature(Map<String, String> protMap, Feature feat) {
            if (feat.getType().contentEquals("CDS")) {
                // Find any useful roles for this feature.
                Collection<Role> featRoles = feat.getUsefulRoles(roleMap);
                if (! featRoles.isEmpty()) {
                    String prot = feat.getProteinTranslation();
                    // Loop through the feature's roles, adding the protein if necessary.
                    for (Role role : featRoles) {
                        String oldRole = protMap.computeIfAbsent(role.getId(), x -> prot);
                        // If this new role is longer, override.  Note that if this is our
                        // first time with the role, the check will fail because of the default.
                        // That is the most common case.
                        if (oldRole.length() < prot.length())
                            protMap.put(role.getId(), prot);
                    }
                }
            }
        }

    }

    @Override
    protected Measurer setupGenome(Genome genome) {
        return this.new Analysis(genome);
    }

    @Override
    protected void parseParms(Map<String, String> keywords) throws ParseFailureException {
        // Get the kmer size.
        this.kSize = this.getIntValue(keywords, "K", 8);
        // Look for a list of roles.
        String roleList = keywords.get("roles");
        if (roleList == null) {
            log.info("Default role list used for protein distance method.");
            this.roleSource = "all";
            this.protRoles = roleMap;
        } else {
            // Get the individual roles and format the method name modifier.
            String[] roles = StringUtils.split(roleList, ',');
            if (roles.length == 1)
                this.roleSource = roles[0];
            else
                this.roleSource = String.format("%s+%d", roles[0], roles.length - 1);
            // Copy the desired roles to the role map.
            this.protRoles = DistanceMethod.customMap(roles);
            log.info("{} roles used for protein distance method.", roles.length);
        }
    }

    @Override
    public double getDistance(Measurer measurer, Measurer other) {
        ProteinDistanceMethod.Analysis m1 = (ProteinDistanceMethod.Analysis) measurer;
        ProteinDistanceMethod.Analysis m2 = (ProteinDistanceMethod.Analysis) other;
        // Process all the roles, accumulating the distances.
        double retVal = this.protRoles.keySet().parallelStream()
                .mapToDouble(x -> this.roleDistance(m1, m2, x)).sum();
        // Divide by the number of roles.
        if (retVal > 0.0)
            retVal /= this.protRoles.size();
        return retVal;
    }

    /**
     * Compute the distance between the two roles for the genomes being measured.  Note
     * that if the role is not common between the genomes, the distance is automatically
     * 1.0 (the maximum).
     *
     * @param m1	analysis of the first genome
     * @param m2	analysis of the second genome
     * @param role	ID of the target role
     *
     * @return the distance between the genomes for that role
     */
    private double roleDistance(Analysis m1, Analysis m2, String role) {
        ProteinKmers k1 = m1.kmerMap.get(role);
        ProteinKmers k2 = m2.kmerMap.get(role);
        double retVal;
        if (k1 == null) {
            if (k2 != null)
                retVal = 1.0;
            else
                retVal = 0.0;
        } else if (k2 == null)
            retVal = 1.0;
        else
            retVal = k1.distance(k2);
        return retVal;
    }

    @Override
    public Details getDetailLevel() {
        return P3Genome.Details.PROTEINS;
    }

    @Override
    public String getName() {
        return String.format("Proteins_%d:%s", this.kSize, this.roleSource);
    }

}
