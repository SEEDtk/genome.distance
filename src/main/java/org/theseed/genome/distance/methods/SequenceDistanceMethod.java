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
import org.theseed.basic.ParseFailureException;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.proteins.Role;
import org.theseed.proteins.RoleMap;
import org.theseed.sequence.ProteinKmers;

/**
 * This is the base class for distance methods that do a kmer comparison of identical roles in genomes.
 * For each role, the longest sequence is kept and converted into kmers.  This gives us greater precision than
 * ProfileDistanceMethod but removes the influence of non-coding regions on DnaDistanceMethod.  Note that the use
 * of the longest sequence simplifies the algorithm, but it leaves us vulnerable to contamination:  a long sequence
 * from an alien species will appear very different from what is supposed to be there.  The subclasses use
 * either DNA sequences or protein sequences.
 *
 * The default for this method is to use all the roles in the standard map; however, a small list of roles can be
 * specified instead.  This allows a single seed protein role (as is used in repgen sets) or a group of common
 * universal roles (as is used in hammer generation).
 *
 * The subclasses will use either DNA sequences or protein sequences, but the ProteinKmers object is used regardless,
 * since the sequences are not reversible.  This is a little jarring, but it is important.
 *
 * The keyword parameters are as follows:
 *
 * K		kmer size
 * roles	comma-delimited list of role IDs (default is to use all)
 * penalty	"Y" to penalize for a role that is present in only one genome, "N" to ignore
 *
 * @author Bruce Parrello
 *
 */
public abstract class SequenceDistanceMethod extends DistanceMethod {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(SequenceDistanceMethod.class);
    /** role map to use */
    private RoleMap protRoles;
    /** protein kmer size */
    protected int kSize;
    /** name for role source */
    protected String roleSource;
    /** TRUE to penalize for a miss */
    protected boolean penalty;

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
            // Loop through the genome features, scanning for useful pegs and storing the sequences.
            genome.getFeatures().parallelStream().forEach(x -> this.scanFeature(protMap, x));
            // Roll up the sequences into a protein kmer map.
            this.kmerMap = new HashMap<String, ProteinKmers>(protMap.size() * 4 / 3);
            for (var protEntry : protMap.entrySet()) {
                ProteinKmers kmers = new ProteinKmers(protEntry.getValue(), SequenceDistanceMethod.this.kSize);
                this.kmerMap.put(protEntry.getKey(), kmers);
            }
        }

        /**
         * Analyze a feature, and if it is a peg with useful roles, store its sequence in the map.
         *
         * @param protMap	target map of role IDs to sequence strings
         * @param feat		feature to analyze
         */
        private void scanFeature(Map<String, String> protMap, Feature feat) {
            if (feat.getType().contentEquals("CDS")) {
                // Find any useful roles for this feature.
                Collection<Role> featRoles = feat.getUsefulRoles(roleMap);
                if (! featRoles.isEmpty()) {
                    String prot = SequenceDistanceMethod.this.getSequence(feat);
                    // Loop through the feature's roles, adding the protein if necessary.
                    for (Role role : featRoles) {
                        String oldRole = protMap.computeIfAbsent(role.getId(), x -> prot);
                        // If this new role is longer, override.  Note that if this is our
                        // first time with the role, the check will fail because the default
                        // is the new role. That is the most common case.
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

    /**
     * @return the sequence to be used for the comparison
     *
     * @param feat	feature whose sequence is desired
     */
    protected abstract String getSequence(Feature feat);

    @Override
    protected void parseParms(Map<String, String> keywords) throws ParseFailureException {
        // Get the kmer size.
        this.kSize = this.getIntValue(keywords, "K", this.getDefaultKmerSize());
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
        // Store the penalty flag.
        this.penalty = this.getBoolValue(keywords, "penalty", this.getPenaltyDefault());
    }


    /**
     * @return TRUE if the default action should be to penalize for roles only in one genome
     */
    protected abstract boolean getPenaltyDefault();

    /**
     * @return the default kmer size
     */
    protected abstract int getDefaultKmerSize();

    @Override
    public double getDistance(Measurer measurer, Measurer other) {
        SequenceDistanceMethod.Analysis m1 = (SequenceDistanceMethod.Analysis) measurer;
        SequenceDistanceMethod.Analysis m2 = (SequenceDistanceMethod.Analysis) other;
        // Process all the roles, accumulating the distances.  Note the distance between a missing role and a present
        // role is always 1.
        int count = 0;
        double retVal = 0.0;
        for (String roleId : this.protRoles.keySet()) {
            ProteinKmers k1 = m1.kmerMap.get(roleId);
            ProteinKmers k2 = m2.kmerMap.get(roleId);
            if (k1 == null) {
                if (k2 != null && this.penalty) {
                    retVal += 1.0;
                    count++;
                }
            } else if (k2 == null) {
                if (this.penalty) {
                    retVal += 1.0;
                    count++;
                }
            } else {
                retVal += k1.distance(k2);
                count++;
            }
        }
        // Divide by the number of roles.
        if (retVal > 0.0)
            retVal /= count;
        return retVal;
    }

    @Override
    public void close() throws Exception {
    }

}
