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
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.ParseFailureException;
import org.theseed.io.TabbedLineReader;
import org.theseed.utils.BasePipeProcessor;
import org.theseed.utils.StringPair;

/**
 * This is a complicated and highly specialized command designed to take an existing genome-comparison test file and augment
 * it to provide more examples of pairs which test distances for specific taxonomic separations, in particular different
 * genera in the same family, different species in the same genus, and different genomes in the same species.
 *
 *
 * The standard input contains the old comparison input file, and the new comparison input file is produced on the standard
 * output.  The old file is presumed to have the genome IDs in the first two columns.  The first column should contain the
 * key genomes.  The new genomes selected will all be compared to at least one existing key genome.  In addition, a complete
 * list of all genomes appearing in at least one output comparison will be written to a specified file.
 *
 * The positional parameter is the name of the evaluation sort file (usually "patric.sort.tbl" in the latest P3Eval directory).
 * The sort file is used to find the new genomes and the GTOs to get the taxonomy of the key genomes.  Only good genomes that
 * have a family, a genus, and a species will be considered.
 *
 * The command-line options are as follows:
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -i	input file name (if not STDIN)
 * -o	output file name (if not STDOUT)
 * -n	number of pairs desired at each taxonomic grouping level (default 1000)
 *
 * --gFile		output file for the full list of genome IDs (default "genomes.tbl" in the current directory)
 *
 * @author Bruce Parrello
 *
 */
public class AugmentProcessor extends BasePipeProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(AugmentProcessor.class);
    /** IDs of primary genomes in the original pairings */
    private Set<String> primaries;
    /** set of pairs already in use */
    private Set<StringPair> pairsUsed;
    /** complete set of genome IDs used */
    private Set<String> genomes;
    /** map of good-genome IDs to genome specs */
    private Map<String, GenomeTaxonSpec> genomeMap;
    /** for each taxonomic grouping level, a map of taxon IDs to sorted genome specs */
    private List<Map<Integer, SortedSet<GenomeTaxonSpec>>> taxonMaps;
    /** count of new pairings needed at each level */
    private int newLeft[];


    // COMMAND-LINE OPTIONS

    /** number of new pairings needed at each level */
    @Option(name = "--num", aliases = { "-n" }, metaVar = "200", usage = "number of pairings needed for each level")
    private int needed;

    /** output file for genome IDs */
    @Option(name = "--gFile", aliases = { "--gfile" }, metaVar = "gList.tbl", usage = "name of output file for genome ID list")
    private File gFile;

    /** evaluation sort file */
    @Argument(index = 0, metaVar = "patric.sort.tbl", usage = "sort file from latest evaluation run", required = true)
    private File sortFile;

    @Override
    protected void setPipeDefaults() {
        this.gFile = new File("genomes.tbl");
        this.needed = 1000;
    }

    @Override
    protected void validatePipeInput(TabbedLineReader inputStream) throws IOException {
        if (inputStream.size() < 2)
            throw new IOException("Input file must have at least two columns.");
    }

    @Override
    protected void validatePipeParms() throws IOException, ParseFailureException {
        // Validate the new-pair count.
        if (this.needed <= 0)
            throw new ParseFailureException("Number of new pairings needed must be greater than 0.");
        // Insure we can read the sort file.
        if (! this.sortFile.canRead())
            throw new FileNotFoundException("Sort file " + this.sortFile + " is not found or unreadable.");
        // We are going to build our main data structures from the sort file, which is read here.
        this.genomeMap = GenomeTaxonSpec.readSortFile(this.sortFile);
        // Create the taxon maps for each taxonomy level.
        this.taxonMaps = new ArrayList<Map<Integer, SortedSet<GenomeTaxonSpec>>>(GenomeTaxonSpec.WORK_LEVELS);
        for (int i = 0; i < GenomeTaxonSpec.WORK_LEVELS; i++)
            this.taxonMaps.add(new HashMap<Integer, SortedSet<GenomeTaxonSpec>>(500));
        // Loop through the specs, storing them.
        for (GenomeTaxonSpec spec : this.genomeMap.values()) {
            // Store this spec in the taxon maps.
            for (int i = 0; i < GenomeTaxonSpec.WORK_LEVELS; i++) {
                var specSet = this.taxonMaps.get(i).computeIfAbsent(spec.getTaxId(i), x -> new TreeSet<GenomeTaxonSpec>());
                specSet.add(spec);
            }
        }
    }

    @Override
    protected void runPipeline(TabbedLineReader inputStream, PrintWriter writer) throws Exception {
        // First we write the headers for the output file.
        writer.println("genome1\tgenome2");
        // Initialize the genome sets.  Here we also count down the pairings we already have at each level.
        this.primaries = new HashSet<String>(3000);
        this.pairsUsed = new HashSet<StringPair>(3000);
        this.genomes = new TreeSet<String>();
        // Clear the counters.  These tell us the number of pairings remaining at each level.  As soon as a count hits
        // zero, we stop searching at that level.
        this.newLeft = new int[GenomeTaxonSpec.WORK_LEVELS];
        Arrays.fill(this.newLeft, this.needed);
        // Now we read in the input file and create the list of primaries and secondaries.  We also count down the
        // distance levels already available.
        log.info("Scanning input file.");
        for (var line : inputStream) {
            String g1 = line.get(0);
            String g2 = line.get(1);
            // We only keep the primaries that have a full taxonomy.
            GenomeTaxonSpec spec1 = this.genomeMap.get(g1);
            if (spec1 != null) {
                this.primaries.add(g1);
                this.pairsUsed.add(new StringPair(g1, g2));
                this.genomes.add(g1);
                this.genomes.add(g2);
                // If the secondary has a full taxonomy, record the distance level.
                GenomeTaxonSpec spec2 = this.genomeMap.get(g2);
                if (spec2 != null) {
                    int lvl = spec1.levelWith(spec2);
                    if (lvl >= 0)
                        this.newLeft[lvl]--;
                }
            }
            // Echo this pair to the output.
            writer.println(g1 + "\t" + g2);
        }
        log.info("{} primary genomes, {} pairs used, {} genomes total.", this.primaries.size(),
                this.pairsUsed.size(), this.genomes.size());
        for (int i = 0; i < GenomeTaxonSpec.WORK_LEVELS; i++) {
            if (newLeft[i] > 0)
                log.info("{} pairs needed for {} level.", newLeft[i], GenomeTaxonSpec.getLevelName(i));
        }
        // Now the main loop.  We loop until we have no needs.
        int found = 0;
        int passes = 0;
        while (this.moreNeeded()) {
            for (var primaryId : this.primaries) {
                // Get the genome spec for this primary.
                GenomeTaxonSpec spec = this.genomeMap.get(primaryId);
                // Try to find new pairings for this genome.
                for (int lvl = 0; lvl < GenomeTaxonSpec.WORK_LEVELS; lvl++) {
                    if (this.newLeft[lvl] > 0) {
                        // Here we need more pairings at this level.  Get the set of other genomes that match.
                        int taxon = spec.getTaxId(lvl);
                        var otherSet = this.taxonMaps.get(lvl).get(taxon);
                        // Only proceed if there are genomes that qualify.  (If there are none, the above returns NULL.)
                        if (otherSet != null) {
                            var iter = otherSet.iterator();
                            // We create a pairing for each genome at the proper taxonomic distance.
                            String newSecondaryId = null;
                            while (newSecondaryId == null && iter.hasNext()) {
                                GenomeTaxonSpec other = iter.next();
                                // Only use this genome if it is not a secondary, not this genome, and it qualifies.
                                String otherId = other.getGenomeId();
                                StringPair testPair = new StringPair(primaryId, otherId);
                                if (! otherId.contentEquals(primaryId) && ! this.pairsUsed.contains(testPair) &&
                                        other.isGoodPairing(spec, lvl)) {
                                    newSecondaryId = otherId;
                                    this.pairsUsed.add(testPair);
                                }
                            }
                            if (newSecondaryId != null) {
                                // Here we found a good pairing. Write the pairing to the main file.
                                writer.println(primaryId + "\t" + newSecondaryId);
                                // Add the other ID to the full genome list.
                                this.genomes.add(newSecondaryId);
                                // Denote one less pairing is needed at this level.
                                this.newLeft[lvl]--;
                                found++;
                            }
                        }
                    }
                }
            }
            passes++;
            log.info("{} passes complete.", passes);
        }
        log.info("{} new pairings found.", found);
        // We have found all the pairings we can.  Update the genome list file.
        try (PrintWriter genomeStream = new PrintWriter(this.gFile)) {
            genomeStream.println("genome_id");
            for (String genomeId : this.genomes)
                genomeStream.println(genomeId);
        }
        log.info("{} genomes required in GTO cache.", this.genomes.size());
    }

    /**
     * @return TRUE if more pairs are needed
     */
    public boolean moreNeeded() {
        return Arrays.stream(this.newLeft).anyMatch(x -> x > 0);
    }

}
