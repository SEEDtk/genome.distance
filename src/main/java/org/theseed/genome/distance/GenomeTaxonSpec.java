/**
 *
 */
package org.theseed.genome.distance;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.io.TabbedLineReader;

/**
 * This object contains a family, genus, species, and genome index arranged in an array.  It is sorted
 * by a floating-point quality score.  We use these objects to find suitable genomes for pairing.
 */
public class GenomeTaxonSpec implements Comparable<GenomeTaxonSpec> {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(GenomeTaxonSpec.class);
    /** quality score */
    private double score;
    /** taxonomic grouping array */
    private final int[] taxons = new int[GenomeTaxonSpec.ALL_LEVELS];
    /** genome ID */
    private String genomeId;
    /** total number of levels */
    public static final int ALL_LEVELS = 4;
    /** number of levels at which we need to pairings */
    public static final int WORK_LEVELS = 3;
    /** GenomeSpec array index for the genome index */
    public static final int GENOME = 3;
    /** GenomeSpec array index for the species ID */
    public static final int SPECIES = 2;
    /** GenomeSpec array index for the genus ID */
    public static final int GENUS = 1;
    /** GenomeSpec array index for the family ID */
    public static final int FAMILY = 0;
    /** list of level names */
    private static final String[] LEVEL_NAMES = new String[] { "family", "genus", "species", "genome" };
    /** next genome index to assign */
    private static int gIndex = 0;

    /**
     * Read the evaluation sort file to create the a genome taxon specification map.
     *
     * @param evalSortFile	evaluation sort file to read
     *
     * @return a map from genome IDs to genome taxon specifications
     *
     * @throws IOException
     */
    public static Map<String, GenomeTaxonSpec> readSortFile(File evalSortFile) throws IOException {
        Map<String, GenomeTaxonSpec> retVal = new HashMap<String, GenomeTaxonSpec>(5000);
        try (TabbedLineReader sortStream = new TabbedLineReader(evalSortFile)) {
            // Find all the columns we need.
            int genomeColIdx = sortStream.findField("Genome");
            int scoreColIdx = sortStream.findField("Score");
            int goodColIdx = sortStream.findField("Good");
            int familyColIdx = sortStream.findField("family");
            int genusColIdx = sortStream.findField("genus");
            int speciesColIdx = sortStream.findField("species");
            log.info("Reading sort file {}.", evalSortFile);
            int readCount = 0;
            for (var line : sortStream) {
                if (line.getFlag(goodColIdx)) {
                    int family = line.getInt(familyColIdx);
                    int genus = line.getInt(genusColIdx);
                    int species = line.getInt(speciesColIdx);
                    // Only process this genome if it has good values for all three taxonomy levels.
                    if (family > 0 && genus > 0 && species > 0) {
                        String genomeId = line.get(genomeColIdx);
                        // Create the genome specification.
                        var spec = new GenomeTaxonSpec(genomeId, line.getDouble(scoreColIdx), family, genus, species);
                        // Store it in the main map.
                        retVal.put(genomeId, spec);
                    }
                }
                readCount++;
                if (log.isInfoEnabled() && readCount % 5000 == 0)
                    log.info("{} genomes processed.", readCount);
            }
            log.info("{} genomes put in master lists.", retVal.size());
        }
        return retVal;
    }


    /**
     * Construct a genome spec for a particular genome.
     *
     * @param genomeId	ID of the genome
     * @param score		quality score
     * @param family	family ID
     * @param genus		genus ID
     * @param species	species ID
     */
    public GenomeTaxonSpec(String genomeId, double score, int family, int genus, int species) {
        this.score = score;
        this.genomeId = genomeId;
        this.taxons[GenomeTaxonSpec.FAMILY] = family;
        this.taxons[GenomeTaxonSpec.GENUS] = genus;
        this.taxons[GenomeTaxonSpec.SPECIES] = species;
        gIndex++;
        this.taxons[GenomeTaxonSpec.GENOME] = gIndex;
    }

    @Override
    public int compareTo(GenomeTaxonSpec o) {
        // Sort first by highest to lowest score.
        int retVal = Double.compare(o.score, this.score);
        // Sort next by genome index (essentially input order from the sort file).
        if (retVal == 0)
            retVal = this.taxons[GenomeTaxonSpec.GENOME] - o.taxons[GenomeTaxonSpec.GENOME];
        return retVal;
    }

    /**
     * @return the quality score
     */
    public double getScore() {
        return this.score;
    }

    /**
     * @return the genomeId
     */
    public String getGenomeId() {
        return this.genomeId;
    }

    /**
     * @return the taxon ID of the appropriate type
     *
     * @param idx	index of the desired grouping level
     */
    public int getTaxId(int idx) {
        return this.taxons[idx];
    }

    /**
     * Determine whether this is a good genome for a pair to represent a grouping at the specified index level.
     *
     * @param spec	GenomeSpec of genome to pair with
     * @param idx	index of the desired grouping level
     */
    public boolean isGoodPairing(GenomeTaxonSpec spec, int idx) {
        return (spec.taxons[idx] == this.taxons[idx] && spec.taxons[idx+1] != this.taxons[idx+1]);
    }

    /**
     * This finds the taxonomic difference level between this genome and another
     *
     * @param spec		other genome's GenomeSpec
     *
     * @return the tightest level at which the genomes are the same, or -1 if they are different at all levels
     */
    public int levelWith(GenomeTaxonSpec spec) {
        int retVal = 0;
        while (retVal < ALL_LEVELS && spec.taxons[retVal] == this.taxons[retVal])
            retVal++;
        return (retVal - 1);
    }

    /**
     * @return the name of a taxonomic difference level
     *
     * @param lvl	level index
     */
    public static String getLevelName(int lvl) {
        String retVal;
        if (lvl < 0)
            retVal = "distant";
        else
            retVal = LEVEL_NAMES[lvl];
        return retVal;
    }

    @Override
    public int hashCode() {
        final int prime = 31;
        int result = 1;
        result = prime * result + ((this.genomeId == null) ? 0 : this.genomeId.hashCode());
        return result;
    }

    @Override
    public boolean equals(Object obj) {
        if (this == obj) {
            return true;
        }
        if (!(obj instanceof GenomeTaxonSpec)) {
            return false;
        }
        GenomeTaxonSpec other = (GenomeTaxonSpec) obj;
        if (this.genomeId == null) {
            if (other.genomeId != null) {
                return false;
            }
        } else if (!this.genomeId.equals(other.genomeId)) {
            return false;
        }
        return true;
    }

}
