/**
 *
 */
package org.theseed.genome.signatures;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.utils.ParseFailureException;

/**
 * This is the base class for all protein classification systems used in creating genome protein signatures.  Its
 * basic function is to take as input a feature and output a set of class names.
 *
 * @author Bruce Parrello
 *
 */
public abstract class SignatureClass {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(SignatureClass.class);
    /** parameter processor */
    protected IParms processor;
    /** expected number of classes per genome */
    protected static final int EXPECTED_CLASSES = 1000;

    /**
     * Interface for controlling classes, used to get additional parameters.
     */
    public interface IParms {

        /**
         * @return the role definition file
         */
        public File getRoleFile();

    }

    /**
     * Enum for signature types
     */
    public static enum Type {
        PGFAM {
            @Override
            public SignatureClass create(IParms processor) {
                return new PgfamSignatureClass(processor);
            }
        }, ROLE {
            @Override
            public SignatureClass create(IParms processor) throws IOException, ParseFailureException {
                return new RoleSignatureClass(processor);
            }
        };

        /**
         * @return a signature classifier of this type
         *
         * @param processor		controlling command processor
         *
         * @throws IOException
         * @throws ParseFailureException
         */
        public abstract SignatureClass create(IParms processor) throws IOException, ParseFailureException;
    }

    /**
     * Construct a signature classifier.
     *
     * @param processor		controlling command processor
     */
    public SignatureClass(IParms processor) {
        this.processor = processor;
    }

    /**
     * Add a feature's classes to the specified set.
     *
     * @param classes	set in which to store the classes
     * @param feat		feature to classify
     */
    protected abstract void addClasses(Set<String> classes, Feature feat);

    /**
     * @return the classes present in a genome
     *
     * @param genome	genome of interest
     */
    public Set<String> getClasses(Genome genome) {
        Set<String> retVal = new HashSet<String>(EXPECTED_CLASSES * 2);
        // Loop through the proteins, collecting the classes found.
        for (Feature feat : genome.getPegs())
            this.addClasses(retVal, feat);
        return retVal;
    }

    /**
     * @return a map of signature class IDs to names
     *
     * @param signatures	collection of signature class IDs to put in the map
     */
    public abstract Map<String, String> getNames(Collection<String> signatures);

}
