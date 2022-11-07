package org.theseed.genome.distance;

import org.junit.jupiter.api.Test;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;

import java.io.File;
import java.io.IOException;

import org.theseed.genome.Genome;
import org.theseed.genome.GenomeDirectory;
import org.theseed.genome.distance.methods.ProtMeasurer;
import org.theseed.proteins.RoleMap;
import org.theseed.proteins.RoleScanner;

/**
 * Unit test for simple App.
 */
public class DistanceTest
{
    /**
     * genome compare test
     *
     * @throws IOException
     * @throws NumberFormatException
     */
    @Test
    public void testComparison() throws NumberFormatException, IOException {
        File rFile = new File("data", "roles.tbl");
        RoleMap roleMap = RoleMap.load(rFile);
        ProtMeasurer.setRoleMap(roleMap);
        File gFile = new File("data", "1121943.4.gto");
        Genome baseGenome = new Genome(gFile);
        ProtMeasurer baseKmers = new ProtMeasurer(baseGenome);
        assertThat(baseKmers.computePercentSimilarity(baseGenome), equalTo(100.0));
        gFile = new File("data", "453962.4.gto");
        Genome otherGenome = new Genome(gFile);
        double closest = baseKmers.computePercentSimilarity(otherGenome);
        assertThat(closest, lessThan(100.0));
        gFile = new File("data", "442341.8.gto");
        otherGenome = new Genome(gFile);
        double closest2 = baseKmers.computePercentSimilarity(otherGenome);
        assertThat(closest2, lessThan(closest));
        gFile = new File("data", "1166948.3.gto");
        otherGenome = new Genome(gFile);
        double closest3 = baseKmers.computePercentSimilarity(otherGenome);
        assertThat(closest3, lessThan(closest2));
        ProtMeasurer otherKmers = new ProtMeasurer(otherGenome);
        double closest3A = otherKmers.computePercentSimilarity(baseGenome);
        assertThat(closest3A, closeTo(closest3, 0.01));
    }

    /**
     * role scanner test
     *
     * @throws IOException
     */
    @Test
    public void testRoleScanner() throws IOException {
        RoleScanner newMap = new RoleScanner();
        assertThat(newMap.fullSize(), equalTo(0));
        File gtoDir = new File("data", "Vipr");
        GenomeDirectory gDir = new GenomeDirectory(gtoDir);
        newMap.addGenomes(gDir);
        assertThat(newMap.fullSize(), equalTo(11));
        assertThat(newMap.getByName("nucleocapsid protein"), not(nullValue()));
        assertThat(newMap.getByName("orf1ab polyprotein"), not(nullValue()));
        assertThat(newMap.getByName("membrane glycoprotein"), not(nullValue()));
    }

}
