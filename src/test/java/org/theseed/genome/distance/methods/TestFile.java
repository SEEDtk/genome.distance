/**
 *
 */
package org.theseed.genome.distance.methods;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;

import java.io.File;
import java.io.IOException;

import org.junit.jupiter.api.Test;
import org.theseed.basic.ParseFailureException;
import org.theseed.genome.Genome;

/**
 * @author Bruce Parrello
 *
 */
class TestFile {

    @Test
    void test() throws ParseFailureException, IOException {
        DistanceMethod dist = DistanceMethod.create("FILE");
        dist.parseParmString("file=data/distances.tbl col1=id1 col2=id2 colx=Proteins_k8:all m0=0 m1=1");
        Genome g1 = new Genome(new File("data", "1007096.3.gto"));
        Genome g2 = new Genome(new File("data", "693746.9.gto"));
        Genome g3 = new Genome(new File("data", "1673721.90.gto"));
        Genome g4 = new Genome(new File("data", "344747.3.gto"));
        Genome g5 = new Genome(new File("data", "2608982.3.gto"));
        Measurer g1m = dist.getMeasurer(g1);
        assertThat(dist.getDistance(g1m, g2), closeTo(0.7493, 0.00005));
        assertThat(dist.getDistance(g1m, g3), closeTo(0.9677, 0.00005));
        Measurer g4m = dist.getMeasurer(g4);
        assertThat(dist.getDistance(g4m, g5), closeTo(0.7537, 0.00005));
    }

}
