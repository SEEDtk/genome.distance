/**
 *
 */
package org.theseed.genome.distance.methods;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;
import java.util.TreeMap;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.ParseFailureException;
import org.theseed.genome.Genome;
import org.theseed.io.TabbedLineReader;
import org.theseed.p3api.P3Genome;

/**
 * This method reads distances from an external method written into a file.  The file must be tab-delimited, with headers
 * The output of the external method needs to be converted into a distance.  Most methods return a fraction or percent
 * similarity.  To convert this to a distance, the user specifies the method value for a distance of 0 and the method
 * value for the distance of 1.  A linear conversion is performed (y = (x - m0) / (m1 - m0)).  The default case is that
 * the method result is a percent similarity, which is the common case.
 *
 * The keywords for this method are as follows.
 *
 * file		the file name
 * name		the name to display for the method (default to the label of the result column)
 * col1		index (1-based) or name of the input column with the first genome ID (default "1")
 * col2		index (1-based) or name of the input column with the second genome ID (default "2")
 * colx		index (1-based) or name of the input column with the method result (default "3")
 * m0		method result value for 0 distance (default 100.0)
 * m1		method result value for 1 distance (default 0.0)
 *
 * @author Bruce Parrello
 *
 */
public class FileDistanceMethod extends DistanceMethod {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(FileDistanceMethod.class);
    /** map of genome IDs to distances maps loaded from the file */
    Map<String, Map<String,Double>> distanceMap;
    /** name for this method */
    private String name;
    /** empty map for nonexistent genomes */
    private static final Map<String, Double> EMPTY_MAP = Collections.emptyMap();

    /**
     * This object contains the distance data for a single genome.  This is a trivial process, since it is all in
     * the map.
     *
     * @param genome	genome to parse
     */
    protected class Analysis extends Measurer {

        /** distance map for this genome */
        private Map<String, Double> gMap;
        /** ID of this genome */
        private String genomeId;

        /**
         * Get the distance map for the specified genome.
         *
         * @param genome	genome whose map is to be loaded
         */
        public Analysis(Genome genome) {
            super(genome);
            this.genomeId = genome.getId();
            this.gMap = FileDistanceMethod.this.distanceMap.getOrDefault(genomeId, EMPTY_MAP);
        }

    }


    @Override
    protected Measurer setupGenome(Genome genome) {
        return this.new Analysis(genome);
    }

    @Override
    protected void parseParms(Map<String, String> keywords) throws ParseFailureException, IOException {
        // This method has to not only parse the keywords but read the file and fill the hash.
        String fileName = keywords.get("file");
        if (fileName == null)
            throw new ParseFailureException("File name required via \"file\" keyword for FILE method.");
        File inFile = new File(fileName);
        if (! inFile.canRead())
            throw new FileNotFoundException("Distance file {} not found or unreadable.");
        // Get the column specs.
        String col1Name = keywords.getOrDefault("col1", "1");
        String col2Name = keywords.getOrDefault("col2", "2");
        String colXName = keywords.getOrDefault("colx", "3");
        // Get the scale factors.
        double m0 = this.getDoubleValue(keywords, "m0", 100.0);
        double m1 = this.getDoubleValue(keywords, "m1", 0.0);
        double scaleFactor = m1 - m0;
        if (Math.abs(scaleFactor) < 1e-10)
            throw new ParseFailureException("m0 and m1 are too close together in FILE method for file " + inFile + ".");
        // Create the master hash.
        this.distanceMap = new HashMap<String, Map<String, Double>>(100);
        // Now read in the file.
        try (TabbedLineReader inStream = new TabbedLineReader(inFile)) {
            int col1Idx = inStream.findField(col1Name);
            int col2Idx = inStream.findField(col2Name);
            int colXIdx = inStream.findField(colXName);
            // Compute the method name here.
            this.name = keywords.getOrDefault("name", inStream.getLabels()[colXIdx]);
            // Loop through the data.
            int count = 0;
            for (var line : inStream) {
                String g1 = line.get(col1Idx);
                String g2 = line.get(col2Idx);
                // Scale the result to a distance.
                double inputVal = line.getDouble(colXIdx);
                double result = (inputVal - m0) / scaleFactor;
                // Store this result in both maps.
                this.storeResult(g1, g2, result);
                this.storeResult(g2, g1, result);
                count++;
            }
            log.info("{} data points read from {}.", count, inFile);
        }
    }

    /**
     * Store a result for the genome pair (g1, g2) in the distance map.
     *
     * @param g1		first genome ID
     * @param g2		second genome ID
     * @param result	distance to store
     */
    private void storeResult(String g1, String g2, double result) {
        var g1Map = this.distanceMap.computeIfAbsent(g1, k -> new TreeMap<String, Double>());
        g1Map.put(g2, result);
    }

    @Override
    public double getDistance(Measurer measurer, Measurer other) {
        FileDistanceMethod.Analysis m1 = (FileDistanceMethod.Analysis) measurer;
        FileDistanceMethod.Analysis m2 = (FileDistanceMethod.Analysis) other;
        // A missing distance counts as 1 (the worst), but we issue a warning.
        Double retVal =  m1.gMap.get(m2.genomeId);
        if (retVal == null) {
            log.warn("WARNING: no {} distance from {} to {}.", this, m1.genomeId, m2.genomeId);
            retVal = 1.0;
        }
        return retVal;
    }

    @Override
    public P3Genome.Details getDetailLevel() {
        return P3Genome.Details.STRUCTURE_ONLY;
    }

    @Override
    public String getName() {
        return this.name;
    }

    @Override
    public void close() throws Exception {
    }

}
