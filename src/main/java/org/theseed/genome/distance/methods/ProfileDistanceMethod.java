/**
 *
 */
package org.theseed.genome.distance.methods;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.IntStream;

import org.theseed.counters.CountMap;
import org.theseed.genome.Genome;
import org.theseed.p3api.P3Genome;
import org.theseed.p3api.P3Genome.Details;
import org.theseed.proteins.RoleMap;

/**
 * This compares role profiles.  The role profile is a vector of role counts, and we return a normalized dot product.
 * Genomes with radically different sequences but the same machinery will compare identical.  We default to the master
 * role map but allow the user to specify an alternate role definition file in the parameters.
 *
 * The parameter string keywords are as follows:
 *
 * roles	file name of alternate role definitions (the file name CANNOT contain spaces or an equal sign)
 *
 * @author Bruce Parrello
 *
 */
public class ProfileDistanceMethod extends DistanceMethod {

    // FIELDS
    /** role map to use */
    private RoleMap profileRoles;
    /** role ID vector used to fill count vectors */
    private List<String> sortedRoles;
    /** vector size */
    private int roleCount;
    /** role file name */
    private String roleSource;

    /**
     * This class contains the data structures for a single genome.
     */
    protected class Analysis extends Measurer {

        /** vector of normalized role counts */
        protected double[] profile;

        /**
         * Analyze a genome to compute its role profile.
         *
         * @param genome	genome to analyze
         */
        public Analysis(Genome genome) {
            super(genome);
            // Count all the role IDs.  CountMap is horribly not thread-safe, so this can't be parallelized.
            CountMap<String> roleCounts = new CountMap<String>();
            genome.getFeatures().stream().filter(x -> x.getType().contentEquals("CDS"))
                    .flatMap(x -> x.getUsefulRoles(ProfileDistanceMethod.this.profileRoles).stream())
                    .forEach(x -> roleCounts.count(x.getId()));
            // Create the output vector.  It is initially all zeroes.
            this.profile = new double[ProfileDistanceMethod.this.roleCount];
            // Now we form the role counts into a vector.  We get the sum of squares to help us normalize.
            double sqSum = roleCounts.counts().stream().mapToDouble(x -> x.getCount() * x.getCount()).sum();
            if (sqSum > 0) {
                // Here the profile is nonzero length, so we can copy it.
                double vecLen = Math.sqrt(sqSum);
                IntStream.range(0,  ProfileDistanceMethod.this.roleCount).parallel()
                        .forEach(i -> this.profile[i] = roleCounts.getCount(ProfileDistanceMethod.this.sortedRoles.get(i)) / vecLen);
            }
        }

    }

    /**
     * This class contains the role profile for a genome.
     */
    @Override
    protected Measurer setupGenome(Genome genome) {
        return this.new Analysis(genome);
    }

    @Override
    protected void parseParms(Map<String, String> keywords) throws IOException {
        // We need to load the role file.  The role source is set here, also.
        String roleFileName = keywords.get("roles");
        if (roleFileName == null) {
            this.profileRoles = DistanceMethod.roleMap;
            this.roleSource = "default";
        } else {
            File roleFile = new File(roleFileName);
            if (! roleFile.canRead())
                throw new FileNotFoundException("Role file " + roleFile + " is not found or unreadable.");
            this.roleSource = roleFile.getName();
            log.info("Loading alternate role file {}.", roleFile);
            this.profileRoles = RoleMap.load(roleFile);
        }
        // Sort the roles to form the scan list.  We use the scan list to fill each vector.
        this.sortedRoles = new ArrayList<String>(this.profileRoles.keySet());
        Collections.sort(this.sortedRoles);
        this.roleCount = this.sortedRoles.size();
        log.info("{} roles in the profile vector.", this.roleCount);
    }

    @Override
    public double getDistance(Measurer measurer, Measurer other) {
        double[] v1 = ((ProfileDistanceMethod.Analysis) measurer).profile;
        double[] v2 = ((ProfileDistanceMethod.Analysis) other).profile;
        // Dot product the two vectors.
        double retVal = 1.0 - IntStream.range(0, this.roleCount).mapToDouble(i -> v1[i] * v2[i]).sum();
        return retVal;
    }

    @Override
    public Details getDetailLevel() {
        return P3Genome.Details.STRUCTURE_ONLY;
    }

    @Override
    public String getName() {
        return String.format("Profile_%s", this.roleSource);
    }

    @Override
    public void close() throws Exception {
    }

}
