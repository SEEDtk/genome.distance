/**
 *
 */
package org.theseed.genome.distance.methods;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.time.Duration;
import java.util.Map;
import java.util.TreeMap;

import org.apache.commons.lang3.StringUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.ParseFailureException;
import org.theseed.genome.Genome;
import org.theseed.p3api.P3Genome;
import org.theseed.proteins.RoleMap;

/**
 * This object represents a measurement method for comparing genomes.
 *
 * Each method must provide a method for parsing additional parameters and for creating an object
 * containing the data required to measure a single genome.  This object will usually be nested so
 * that it has easy access to the parent method.  A master role map is maintained, and
 * utility methods provided for extracting a sub-list of useful roles.
 *
 * @author Bruce Parrello
 *
 */
public abstract class DistanceMethod implements AutoCloseable {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(DistanceMethod.class);
    /** role definition map */
    protected static RoleMap roleMap;
    /** name of this method */
    private String methodName;

    /**
     * This enum represents the method type.
     */
    public enum Type {
        /** protein-role comparison */
        PROTEIN {
            @Override
            public DistanceMethod create() {
                return new ProteinDistanceMethod();
            }
        },
        /** dna-role comparison */
        GENE {
            @Override
            public DistanceMethod create() {
                return new GeneDistanceMethod();
            }
        },
        /** role profile comparison */
        PROFILE {
            @Override
            public DistanceMethod create() {
                return new ProfileDistanceMethod();
            }
        },
        /** whole-genome DNA kmer distance */
        DNA {
            @Override
            public DistanceMethod create() {
                return new DnaDistanceMethod();
            }
        },
        /** average nucleotide identity */
        ANI {
            @Override
            public DistanceMethod create() {
                return new AniDistanceMethod();
            }
        },
        /** SSU kmer distance */
        SSU {
            @Override
            public DistanceMethod create() {
                return new SsuDistanceMethod();
            }
        },
        /** crude distance based on first differing taxonomic group */
        TAXONOMY {
            @Override
            public DistanceMethod create() {
                return new TaxonDistanceMethod();
            }
        },
        /** distance loaded from a file of pre-computed values */
        FILE {
            @Override
            public DistanceMethod create() {
                return new FileDistanceMethod();
            }
        };

        /**
         * @return a method descriptor of this type
         */
        public abstract DistanceMethod create();
    }

    /**
     * Compute a distance method from a type name.
     *
     * @param type	method type name
     *
     * @throws ParseFailureException
     */
    public static DistanceMethod create(String type) throws ParseFailureException {
        DistanceMethod retVal;
        try {
            Type t = Type.valueOf(type.toUpperCase());
            retVal = t.create();
        } catch (IllegalArgumentException e) {
            throw new ParseFailureException("Invalid distance method type " + type + ".");
        }
        return retVal;
    }

    /**
     * Load the master role map.
     *
     * @param roleFile	role definition file
     *
     * @throws IOException
     */
    public static void loadRoles(File roleFile) throws IOException {
        if (! roleFile.canRead())
            throw new FileNotFoundException("Role file " + roleFile + " is not found or unreadable.");
        roleMap = RoleMap.load(roleFile);
        log.info("{} roles from {} loaded into master role map.", roleMap.size(), roleFile);
    }

    /**
     * Create a measurer that uses this method for a specified genome.  This public method
     * emits useful tracing information.
     *
     * @param genome	genome to parse into measurement data
     *
     * @return a measurer for a specified genome
     */
    public Measurer getMeasurer(Genome genome) {
        long start = System.currentTimeMillis();
        log.debug("Processing method {} for genome {}.", this, genome);
        Measurer retVal = this.setupGenome(genome);
        if (log.isInfoEnabled()) {
            var duration = Duration.ofMillis(System.currentTimeMillis() - start);
            log.debug("{} to process genome of length {}.", duration, genome.getLength());
        }
        return retVal;
    }

    /**
     * Create a measurer for a specific genome that uses this method.
     *
     * @param genome	genome to parse into measurement data
     *
     * @return the created measurer
     */
    protected abstract Measurer setupGenome(Genome genome);

    /**
     * Parse the parameter string for this method to extract its tuning parameters.  The parameters are
     * in the form of key-value pairs consisting of a keyword, an equal sign, and a value.  The pairs
     * are delimited by whitespace.
     *
     * @param parse		parameter string to parse
     *
     * @throws IOException
     * @throws ParseFailureException
     */
    public void parseParmString(String parse) throws IOException, ParseFailureException {
        // Build the keyword map.
        Map<String, String> keywords = new TreeMap<String, String>();
        String[] parms = StringUtils.split(parse);
        for (String parm : parms) {
            int sep = parm.indexOf('=');
            if (sep <= 0)
                throw new ParseFailureException("Invalid keyword string \"" + parm + "\" in measurement method parameters.");
            keywords.put(StringUtils.substring(parm, 0, sep), StringUtils.substring(parm, sep+1));
        }
        // Ask the subclass to process the parameters.
        this.parseParms(keywords);
        // Now save the method name.
        this.methodName = this.getName();
    }

    /**
     * Parse the parameters for this measurement method.
     *
     * @param keywords	key-value pair mapping for the method parameters
     *
     * @throws ParseFailureException
     * @throws IOException
     */
    protected abstract void parseParms(Map<String, String> keywords) throws ParseFailureException, IOException;

    /**
     * Compare two genomes.  In this case the first genome has an existing measurer and the second does
     * not.
     *
     * @param	measurer	measurer using this method for a genome
     * @param	genome		genome to compare
     *
     * @return the distance between the two genomes (from 0 to 1)
     */
    public double getDistance(Measurer measurer, Genome genome) {
        // Get a measurer for the other genome.
        Measurer other = this.getMeasurer(genome);
        // Compute the distance.
        double retVal = this.getDistance(measurer, other);
        return retVal;
    }

    /**
     * This computes the distance using two pre-compiled measurers.  They must be of the
     * appropriate type.
     *
     * @param measurer		measurer for the first genome
     * @param other			measurer for the other genome
     *
     * @return the distance between the genomes (from 0 to 1)
     */
    public abstract double getDistance(Measurer measurer, Measurer other);

    /**
     * Compute a reduced role map containing a subset of the roles in the master map.
     * This is optimized for the most common case, which is a single role; however,
     * since it is only done once per method, the overhead is minimal even for multiple
     * roles.
     *
     * @param roles		list of IDs for roles to keep
     *
     * @return a reduced role map
     *
     * @throws ParseFailureException
     */
    protected static RoleMap customMap(String... roles) throws ParseFailureException {
        RoleMap retVal = new RoleMap();
        for (String role : roles) {
            var roleObjects = roleMap.getAllById(role);
            if (roleObjects.isEmpty())
                throw new ParseFailureException("Role " + role + " not found in master role map.");
            retVal.putAll(roleObjects);
        }
        return retVal;
    }

    /**
     * @return the detail level required for genomes using this method
     */
    public abstract P3Genome.Details getDetailLevel();

    /**
     * @return the abbreviated name of this method (suitable for column headings and tracing)
     */
    public abstract String getName();

    @Override
    public String toString() {
        String retVal;
        if (this.methodName == null)
            retVal = this.getName();
        else
            retVal = this.methodName;
        return retVal;
    }

    /**
     * This is a utility method for processing a boolean (Y/N) keyword.
     *
     * @param keywords			keyword map
     * @param key				keyword name
     * @param defaultVal		default value
     *
     * @return the boolean value
     *
     * @throws ParseFailureException
     */
    public boolean getBoolValue(Map<String, String> keywords, String key, boolean defaultVal) throws ParseFailureException {
        String val = keywords.get(key);
        boolean retVal;
        if (val == null)
            retVal = defaultVal;
        else if (val.contentEquals("Y"))
            retVal = true;
        else if (val.contentEquals("N"))
            retVal = false;
        else
            throw new ParseFailureException("Invalid Y/N flag for keyword " + key + ".");
        return retVal;
    }

    /**
     * This is a utility method for processing an integer keyword.
     *
     * @param keywords		keyword map
     * @param key			keyword name
     * @param defaultVal	default value
     *
     * @return the integer value
     *
     * @throws ParseFailureException
     */
    public int getIntValue(Map<String, String> keywords, String key, int defaultVal) throws ParseFailureException {
        String val = keywords.get(key);
        int retVal;
        if (val == null)
            retVal = defaultVal;
        else try {
            retVal = Integer.valueOf(val);
        } catch (NumberFormatException e) {
            throw new ParseFailureException("Invalid numeric for keyword " + key + ".");
        }
        return retVal;
    }

    /**
     * This is a utility method for processing a floating-point keyword.
     *
     * @param keywords		keyword map
     * @param key			keyword name
     * @param defaultVal	default value
     *
     * @return the integer value
     *
     * @throws ParseFailureException
     */
    public double getDoubleValue(Map<String, String> keywords, String key, double defaultVal) throws ParseFailureException {
        String val = keywords.get(key);
        double retVal;
        if (val == null)
            retVal = defaultVal;
        else try {
            retVal = Double.valueOf(val);
        } catch (NumberFormatException e) {
            throw new ParseFailureException("Invalid numeric for keyword " + key + ".");
        }
        return retVal;
    }

}
