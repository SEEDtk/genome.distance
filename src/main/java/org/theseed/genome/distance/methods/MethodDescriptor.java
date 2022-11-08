/**
 *
 */
package org.theseed.genome.distance.methods;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Map;
import java.util.TreeMap;

import org.apache.commons.lang3.StringUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Genome;
import org.theseed.proteins.RoleMap;
import org.theseed.utils.ParseFailureException;

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
public abstract class MethodDescriptor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(MethodDescriptor.class);
    /** role definition map */
    private static RoleMap roleMap;

    /**
     * This enum represents the method type.
     */
    public enum Type {
        /** seed protein kmer distance */
        SEED {
            @Override
            public MethodDescriptor create() {
                // TODO code for create SEED
                return null;
            }
        },
        /** protein-role comparison */
        PROTEIN {
            @Override
            public MethodDescriptor create() {
                // TODO code for create PROTEIN
                return null;
            }
        },
        /** whole-genome DNA kmer distance */
        DNA {
            @Override
            public MethodDescriptor create() {
                // TODO code for create DNA
                return null;
            }
        },
        /** SSU kmer distance */
        SSU {
            @Override
            public MethodDescriptor create() {
                // TODO code for create SSU
                return null;
            }
        };

        /**
         * @return a method descriptor of this type
         */
        public abstract MethodDescriptor create();
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
     * Create a measurer for a specific genome that uses this method.
     *
     * @param genome	genome to parse into measurement data
     */
    public abstract Measurer getMeasurer(Genome genome);

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
        Map<String, String> keywords = new TreeMap<String, String>();
        String[] parms = StringUtils.split(parse);
        for (String parm : parms) {
            int sep = parm.indexOf('=');
            if (sep <= 0)
                throw new ParseFailureException("Invalid keyword string \"" + parm + "\" in measurement method parameters.");
            keywords.put(StringUtils.substring(parm, 0, sep), StringUtils.substring(parm, sep+1));
        }
        this.parseParms(keywords);
    }

    /**
     * Parse the parameters for this measurement method.
     *
     * @param keywords	key-value pair mapping for the method parameters
     */
    protected abstract void parseParms(Map<String, String> keywords);

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

}
