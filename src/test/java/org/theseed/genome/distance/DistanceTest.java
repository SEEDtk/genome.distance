package org.theseed.genome.distance;

import junit.framework.Test;

import junit.framework.TestCase;
import junit.framework.TestSuite;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;

import java.io.File;
import java.io.IOException;

import org.theseed.genome.Genome;
import org.theseed.genome.GenomeDirectory;
import org.theseed.proteins.RoleMap;
import org.theseed.proteins.RoleScanner;

/**
 * Unit test for simple App.
 */
public class DistanceTest extends TestCase
{
    /**
     * Create the test case
     *
     * @param testName name of the test case
     */
    public DistanceTest( String testName )
    {
        super( testName );
    }

    /**
     * @return the suite of tests being tested
     */
    public static Test suite()
    {
        return new TestSuite( DistanceTest.class );
    }

    /**
     * genome compare test
     *
     * @throws IOException
     * @throws NumberFormatException
     */
    public void testComparison() throws NumberFormatException, IOException {
        File rFile = new File("src/test", "roles.tbl");
        RoleMap roleMap = RoleMap.load(rFile);
        File gFile = new File("src/test", "1166948.3.gto");
        Genome baseGenome = new Genome(gFile);
        Measurer baseKmers = new Measurer(baseGenome, roleMap);
        assertThat(baseKmers.computePercentSimilarity(baseGenome), equalTo(100.0));
    }

    /**
     * role scanner test
     *
     * @throws IOException
     */
    public void testRoleScanner() throws IOException {
        RoleScanner newMap = new RoleScanner();
        assertThat(newMap.fullSize(), equalTo(0));
        File gtoDir = new File("src/test", "Vipr");
        GenomeDirectory gDir = new GenomeDirectory(gtoDir);
        newMap.addGenomes(gDir);
        assertThat(newMap.fullSize(), equalTo(11));
        assertNotNull(newMap.getByName("nucleocapsid protein"));
        assertNotNull(newMap.getByName("orf1ab polyprotein"));
        assertNotNull(newMap.getByName("membrane glycoprotein"));
    }

}
