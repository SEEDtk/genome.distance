/**
 *
 */
package org.theseed.genome.distance.methods;

import java.io.IOException;

import org.theseed.genome.Genome;
import org.theseed.p3api.P3Genome;
import org.theseed.p3api.P3Genome.Details;
import org.theseed.proteins.RoleMap;
import org.theseed.sequence.GenomeKmers;
import org.theseed.utils.ParseFailureException;

/**
 * This is the base class for a genome distance measurer.  The measurer is constructed from
 * the genome object itself.  It contains the genome ID and name at this level.  The subclass
 * must also contain information that can be used to compute the percent similarity between
 * genomes.
 *
 * @author Bruce Parrello
 *
 */
public abstract class Measurer {

    // FIELDS
    /** genome ID */
    private String genomeId;
    /** genome name */
    private String genomeName;

    /**
     * Enumeration for the different measurement types.
     */
    public static enum Type {
        /** kmer distance on a protein/protein basis */
        PROTEIN {
            @Override
            public Measurer create(Genome genome) {
                return new ProtMeasurer(genome);
            }

            @Override
            public void init(IParms processor) throws ParseFailureException {
                // Here we set up the role definitions.
                var roleMap = processor.getRoleMap();
                if (roleMap == null)
                    throw new ParseFailureException("Role definitions required for PROTEIN measurement.");
                ProtMeasurer.setRoleMap(roleMap);
            }

            @Override
            protected P3Genome.Details getLevel() {
                return P3Genome.Details.PROTEINS;
            }

        },

        /** raw DNA kmer distance */
        CONTIG {
            @Override
            public Measurer create(Genome genome) {
                return new DnaMeasurer(genome);
            }

            @Override
            public void init(IParms processor) throws ParseFailureException {
                // Here we set up the DNA kmer size.
                int kSize = processor.getKmerSize();
                if (kSize <= 0)
                    throw new ParseFailureException("Kmer size must be at least 1 for CONTIG measurement.");
                GenomeKmers.setKmerSize(kSize);
            }

            @Override
            protected P3Genome.Details getLevel() {
                return P3Genome.Details.CONTIGS;
            }

        },

        /** kmer distance for a single seed protein */
        SEED {

            @Override
            public Measurer create(Genome genome) {
                return new SeedMeasurer(genome);
            }

            @Override
            public void init(IParms processor) throws IOException, ParseFailureException {
                // Here we set up the role definition.
                var roleMap = processor.getRoleMap();
                if (roleMap == null)
                    throw new ParseFailureException("Role definitions required for SEED measurement.");
                String seedId = processor.getSeedId();
                if (seedId == null)
                    throw new ParseFailureException("Seed protein ID required for SEED measurement.");
                SeedMeasurer.setRole(roleMap, seedId);
            }

            @Override
            protected Details getLevel() {
                return P3Genome.Details.PROTEINS;
            }

        },

        /** kmer distance for SSUs */
        SSU {

            @Override
            public Measurer create(Genome genome) {
                return new SsuMeasurer(genome);
            }

            @Override
            public void init(IParms processor) throws IOException, ParseFailureException {
                // Here we set up the DNA kmer size.
                int kSize = processor.getKmerSize();
                if (kSize <= 0)
                    throw new ParseFailureException("Kmer size must be at least 1 for SSU measurement.");
                GenomeKmers.setKmerSize(kSize);
            }

            @Override
            protected Details getLevel() {
                return P3Genome.Details.FULL;
            }

        }
        ;

        /**
         * @return a measurer of this type for the specified genome.
         *
         * @param genome	genome of interest
         */
        public abstract Measurer create(Genome genome);

        /**
         * Initialize measurement for the specified command processor.
         *
         * @params processor	controlling command processor
         */
        public abstract void init(IParms processor) throws IOException, ParseFailureException;

        /**
         * @return the genome detail level required to perform this comparison
         */
        protected abstract P3Genome.Details getLevel();
    }

    /**
     * Interface for processors using a measurer.
     */
    public interface IParms {

        /**
         * @return the role definition map
         */
        public RoleMap getRoleMap();

        /**
         * @return the DNA kmer size
         */
        public int getKmerSize();

        /**
         * @return the seed protein ID
         */
        public String getSeedId();

    }

    /**
     * Construct a measurer.
     *
     * @param genome	genome to be measured for distances
     */
    public Measurer(Genome genome) {
        this.genomeId = genome.getId();
        this.genomeName = genome.getName();
    }

    /**
     * @return the genome ID
     */
    public String getId() {
        return this.genomeId;
    }

    /**
     * @return the genome name
     */
    public String getName() {
        return this.genomeName;
    }

    /**
     * @return the type of this measurer
     */
    public abstract Measurer.Type getType();

    /**
     * @return the percent similarity between two genomes
     *
     * @param other		other genome to compare to this one
     */
    public double computePercentSimilarity(Genome other) {
        Measurer otherMeasurer = this.getType().create(other);
        return this.computePercentSimilarity(otherMeasurer);
    }

    /**
     * This computes the percent similarity between this measurer and the measurer for another
     * genome.  The other measurer must be of the same type.
     *
     * @param otherMeasurer		measurer for the other genome
     *
     * @return the percent similarity between the two genomes
     */
    public abstract double computePercentSimilarity(Measurer otherMeasurer);

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
        if (!(obj instanceof Measurer)) {
            return false;
        }
        Measurer other = (Measurer) obj;
        if (this.genomeId == null) {
            if (other.genomeId != null) {
                return false;
            }
        } else if (!this.genomeId.equals(other.genomeId)) {
            return false;
        }
        return true;
    }

    /**
     * @return the detail level required for genomes when using the specified comparison
     *
     * @param comparisonType	type of comparison to check
     */
    public static Details getLevel(Type comparisonType) {
        return comparisonType.getLevel();
    }

}
