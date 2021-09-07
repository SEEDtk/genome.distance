/**
 *
 */
package org.theseed.genome.distance;

import java.io.File;
import java.io.IOException;

import org.theseed.genome.Genome;
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
        PROTEIN {
            @Override
            public Measurer create(Genome genome) {
                return new ProtMeasurer(genome);
            }

            @Override
            public void init(IParms processor) throws ParseFailureException {
                // Here we set up the role definitions.
                File roleFile = processor.getRoleFile();
                if (roleFile == null)
                    throw new ParseFailureException("Role definition file required for PROTEIN measurement.");
                RoleMap roles = RoleMap.load(roleFile);
                ProtMeasurer.setRoleMap(roles);
            }

        }, CONTIG {
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

        };

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
    }

    /**
     * Interface for processors using a measurer.
     */
    public interface IParms {

        /**
         * @return the role definition file
         */
        public File getRoleFile();

        /**
         * @return the DNA kmer size
         */
        public int getKmerSize();

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
     * @return the percent similarity between two genomes
     *
     * @param other		other genome to compare to this one
     */
    public abstract double computePercentSimilarity(Genome other);

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

}
