/**
 *
 */
package org.theseed.genome.distance.methods;

import java.io.File;
import java.io.IOException;

import org.junit.jupiter.api.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.ParseFailureException;
import org.theseed.genome.Genome;

/**
 * @author Bruce Parrello
 */
class TestAni {

    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(TestAni.class);


    @Test
    void test() throws IOException, ParseFailureException {
        Genome g1 = new Genome(new File("data", "1007096.3.gto"));
        Genome g2 = new Genome(new File("data", "693746.9.gto"));
        Genome g3 = new Genome(new File("data", "1673721.90.gto"));
        Genome g4 = new Genome(new File("data", "344747.3.gto"));
        Genome g5 = new Genome(new File("data", "2608982.3.gto"));
        Genome g6 = new Genome(new File("data", "1140.7.gto"));
        Genome g7 = new Genome(new File("data", "1357541.3.gto"));
        // Known values.  g1/g2 = 81.72%, g1/g3 = 69.9%, g4/g5 = 74.8%, g6/g7 = 99.98%.
        DistanceMethod ani = DistanceMethod.create("ANI");
        ani.parseParmString("");
        Measurer g1m = ani.getMeasurer(g1);
        double g1_2 = ani.getDistance(g1m, g2);
        double g1_3 = ani.getDistance(g1m, g3);
        Measurer g4m = ani.getMeasurer(g4);
        double g4_5 = ani.getDistance(g4m, g5);
        Measurer g6m = ani.getMeasurer(g6);
        double g6_7 = ani.getDistance(g6m, g7);
        log.info("g1 to g2 is {}, {}%.", g1_2, (1.0 - g1_2) * 100);
        log.info("g1 to g3 is {}, {}%.", g1_3, (1.0 - g1_3) * 100);
        log.info("g4 to g5 is {}, {}%.", g4_5, (1.0 - g4_5) * 100);
        log.info("g6 to g7 is {}, {}%.", g6_7, (1.0 - g6_7) * 100);
    }

}
