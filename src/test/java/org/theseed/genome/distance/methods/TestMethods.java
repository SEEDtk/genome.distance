/**
 *
 */
package org.theseed.genome.distance.methods;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.junit.jupiter.api.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.ParseFailureException;
import org.theseed.genome.Genome;

/**
 * @author Bruce Parrello
 *
 */
public class TestMethods {

    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(TestMethods.class);
    /** test genome IDs, in order of increasing distance */
    public static final String[] TEST_GENOMES = new String[] { "1062.144", "366616.15", "2720031.3", "272844.11" };
    /** main genome ID */
    public static final String BASE_GENOME = "1062.102";
    /** list of method names */
    public static final String[] METHODS = { "PROFILE", "DNA", "PROTEIN", "PROTEIN", "PROTEIN", "SSU", "ANI", "GENE", "TAXONOMY" };
    /** list of parameter strings */
    public static final String[] PARMS = {
            /* PROFILE */	"roles=data/roles.to.use",
            /* DNA */ 		"K=20",
            /* PROTEIN */	"",
            /* PROTEIN */	"roles=PolyNucl,DnaDireRnaPolyBeta3 K=8",
            /* PROTEIN */	"roles=PhenTrnaSyntAlph",
            /* SSU */		"K=20",
            /* ANI */		"",
            /* GENE */		"",
            /* TAXONOMY */	""
        };

    @Test
    void testAllMethods() throws IOException, ParseFailureException {
        File roleFile = new File("data", "roles.in.subsystems");
        DistanceMethod.loadRoles(roleFile);
        Genome baseGenome = new Genome(new File("data", BASE_GENOME + ".gto"));
        List<Genome> testGenomes = new ArrayList<Genome>(TEST_GENOMES.length);
        log.info("Base genome is {}.", baseGenome);
        for (String genomeId : TEST_GENOMES) {
            var genome = new Genome(new File("data", genomeId + ".gto"));
            testGenomes.add(genome);
        }
        // Now we have all the genomes.  Loop through the methods, creating measurers.
        for (int i = 0; i < METHODS.length; i++) {
            DistanceMethod method = DistanceMethod.create(METHODS[i]);
            method.parseParmString(PARMS[i]);
            Measurer baseAnalysis = method.getMeasurer(baseGenome);
            double lastDist = method.getDistance(baseAnalysis, baseAnalysis);
            assertThat(lastDist, closeTo(0.0, 1e-10));
            for (Genome testGenome : testGenomes) {
                double newDist = method.getDistance(baseAnalysis, testGenome);
                assertThat(newDist, lessThanOrEqualTo(1.0));
                assertThat(newDist, greaterThanOrEqualTo(0.0));
                log.info("Distance from base to {} using {} is {}.", testGenome, method, newDist);
                lastDist = newDist;
            }
        }
    }

}
