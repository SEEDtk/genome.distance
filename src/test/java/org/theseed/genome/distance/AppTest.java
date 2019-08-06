package org.theseed.genome.distance;

import junit.framework.Test;

import junit.framework.TestCase;
import junit.framework.TestSuite;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;

import java.io.File;
import java.io.IOException;

import org.theseed.genome.Genome;
import org.theseed.proteins.RoleMap;

/**
 * Unit test for simple App.
 */
public class AppTest extends TestCase
{
    /**
     * Create the test case
     *
     * @param testName name of the test case
     */
    public AppTest( String testName )
    {
        super( testName );
    }

    /**
     * @return the suite of tests being tested
     */
    public static Test suite()
    {
        return new TestSuite( AppTest.class );
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

}
