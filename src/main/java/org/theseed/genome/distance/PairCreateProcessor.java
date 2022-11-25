/**
 *
 */
package org.theseed.genome.distance;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.stream.IntStream;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.io.TabbedLineReader;
import org.theseed.utils.BaseReportProcessor;
import org.theseed.utils.ParseFailureException;

/**
 * This command will create a list of genome pairs to test, attempting to balance as much as possible the
 * spread across the tree and at different distances.  To do this, it reads a repgen list file and organizes
 * the genomes by distance within representative.  For each representative, we attempt to find a neighbor
 * at the minimum distance, the maximum distance, the median distance, and the quarter-median distance.
 *
 * The positional parameters are the name of the repgen list file and the name of an evaluation sort file.
 * The former is the primary input; the latter is used to verify that the genomes are good and that they
 * have a known family, genus, and species.
 *
 * The command-line options are as follows.
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -o	output file for the pairings (if not STDOUT)
 *
 * --gFile	name of a file to contain the list of genome IDs used (default "genomes.tbl" in the current directory)
 *
 * @author Bruce Parrello
 *
 */
public class PairCreateProcessor extends BaseReportProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(PairCreateProcessor.class);
    /** map of genome IDs to taxonomy specification objects */
    private Map<String, GenomeTaxonSpec> genomeMap;
    /** map of representatives to neighbor lists */
    private Map<String, List<Neighbor>> neighborhoodMap;
    /** number of pairs output */
    private int outCount;
    /** set containing the IDs of the genomes used */
    private Set<String> genomes;
    /** number of pairs at each taxonomic distance */
    private int[] taxDistanceCounts;
    /** number of distant genome pairs */
    private int farPairCount;

    // COMMAND-LINE OPTIONS

    /** output file for genome IDs */
    @Option(name = "--gFile", aliases = { "--gfile" }, metaVar = "gList.tbl", usage = "name of output file for genome ID list")
    private File gFile;

    /** name of the repgen list file */
    @Argument(index = 0, metaVar = "repgen.list.tbl", usage = "name of the repgen list file containing the neighbor data", required = true)
    private File repListFile;

    /** name of the evaluation sort file */
    @Argument(index = 1, metaVar = "patric.sort.tbl", usage = "name of the evaluation sort table containing taxonomic specs", required = true)
    private File sortFile;

    /**
     * This object contains a neighbor genome ID, its similarity score, and its distance from the
     * representative. It is sorted from closest genome to furthest.
     */
    protected static class Neighbor implements Comparable<Neighbor> {

        /** ID of the neighbor genome */
        private String genomeId;
        /** distance to the neighbor from the representative */
        private double distance;
        /** similarity score between the neighbor and the representative */
        private int simScore;

        /**
         * Create a neighbor record.
         *
         * @param gId		ID of the neighbor genome
         * @param dist		distance to the representative
         * @param sim		similarity score with the representative
         */
        public Neighbor(String gId, double dist, int sim) {
            this.genomeId = gId;
            this.distance = dist;
            this.simScore = sim;
        }

        /**
         * @return the neighbor's genome ID
         */
        public String getId() {
            return this.genomeId;
        }

        @Override
        public int compareTo(Neighbor o) {
            int retVal = Double.compare(this.distance, o.distance);
            if (retVal == 0) {
                // A higher similarity score means a closer neighbor.
                retVal = o.simScore - this.simScore;
                if (retVal == 0)
                    retVal = this.genomeId.compareTo(o.genomeId);
            }
            return retVal;
        }

    }

    @Override
    protected void setReporterDefaults() {
        this.gFile = new File("genomes.tbl");
    }

    @Override
    protected void validateReporterParms() throws IOException, ParseFailureException {
        // Verify that both input files are readable.
        if (! this.sortFile.canRead())
            throw new FileNotFoundException("Evaluation sort file " + this.sortFile + " is not found or unreadable.");
        if (! this.repListFile.canRead())
            throw new FileNotFoundException("Repgen list file " + this.repListFile + " is not found or unreadable.");
        // Now read in the sort file to get the good genomes.
        this.genomeMap = GenomeTaxonSpec.readSortFile(this.sortFile);
        // Create the genome save set.  We use a sorted set so they come out in a user-friendly order.
        this.genomes = new TreeSet<String>();
    }

    @Override
    protected void runReporter(PrintWriter writer) throws Exception {
        // Here we read the repgen list file to get the neighbor data.  Note that only genomes in the genome map are kept.
        this.neighborhoodMap = new HashMap<String, List<Neighbor>>(2000);
        log.info("Reading repgen list file {}.", this.repListFile);
        try (TabbedLineReader repStream = new TabbedLineReader(this.repListFile)) {
            // Find the important columns.
            int neighborColIdx = repStream.findField("genome_id");
            int repColIdx = repStream.findField("rep_id");
            int distColIdx = repStream.findField("distance");
            int simColIdx = repStream.findField("score");
            // Now loop through the file.
            int count = 0;
            int stored = 0;
            for (var line : repStream) {
                String neighborId = line.get(neighborColIdx);
                String repId = line.get(repColIdx);
                // Only proceed if both genomes are good and they are not equal.
                if (this.genomeMap.containsKey(neighborId) && this.genomeMap.containsKey(repId) &&
                        ! neighborId.contentEquals(repId)) {
                    // Now we build the neighbor object and attach it to the representative's neighbor list.
                    var neighbor = new Neighbor(neighborId, line.getDouble(distColIdx), line.getInt(simColIdx));
                    var neighborList = this.neighborhoodMap.computeIfAbsent(repId, x -> new ArrayList<Neighbor>());
                    neighborList.add(neighbor);
                    stored++;
                }
                count++;
                if (count % 5000 == 0)
                    log.info("{} genomes processed.  {} stored.", count, stored);
            }
        }
        log.info("{} representatives are good and have a neighborhood.", this.neighborhoodMap.size());
        // Now we produce the output.  Start with the header.
        writer.println("genome1\tgenome2");
        // Set up the counters.
        this.taxDistanceCounts = new int[GenomeTaxonSpec.WORK_LEVELS];
        Arrays.fill(this.taxDistanceCounts, 0);
        this.outCount = 0;
        this.farPairCount = 0;
        // Loop through the neighborhoods, producing output.
        for (var neighborhoodEntry : this.neighborhoodMap.entrySet()) {
            String repId = neighborhoodEntry.getKey();
            List<Neighbor> neighborhood = neighborhoodEntry.getValue();
            // If there are four or fewer points, keep them all.
            final int nSize = neighborhood.size();
            if (nSize <= 4) {
                for (var neighbor : neighborhood)
                    this.output(writer, repId, neighbor);
                log.info("{} pairs output for {}.", nSize, repId);
            } else {
                // Sort the neighborhood list.
                Collections.sort(neighborhood);
                // Get the four desired points.  Because the neighborhood is big enough, they will be different.
                 int[] points = new int[] { 0, (nSize / 4), (nSize / 2), nSize - 1 };
                 Arrays.stream(points).forEach(i -> this.output(writer, repId, neighborhood.get(i)));
                 log.info("Standard pairs output for {}.", repId);
            }
        }
        log.info("{} total pairs output.", this.outCount);
        IntStream.range(0, GenomeTaxonSpec.WORK_LEVELS)
            .forEach(i -> log.info("{} pairs at {} level.", this.taxDistanceCounts[i], GenomeTaxonSpec.getLevelName(i)));
        log.info("{} pairs are far apart.", this.farPairCount);
        // Now write the genome ID list.
        log.info("{} genomes required.  Writing list.", this.genomes.size());
        try (PrintWriter gWriter = new PrintWriter(this.gFile)) {
            gWriter.println("genome_id");
            for (var genomeId : this.genomes)
                gWriter.println(genomeId);
        }
    }

    /**
     * Write out a single genome pair.
     *
     * @param writer		output stream
     * @param repId			representative genome ID
     * @param neighbor		neighbor object for the second genome
     */
    private void output(PrintWriter writer, String repId, Neighbor neighbor) {
        final String genomeId = neighbor.getId();
        writer.println(repId + "\t" + genomeId);
        this.outCount++;
        this.genomes.add(repId);
        this.genomes.add(genomeId);
        // Update the taxonomic distance counts.
        var spec1 = this.genomeMap.get(repId);
        var spec2 = this.genomeMap.get(genomeId);
        int lvl = spec1.levelWith(spec2);
        if (lvl >= 0)
            this.taxDistanceCounts[lvl]++;
        else
            this.farPairCount++;
    }

}
