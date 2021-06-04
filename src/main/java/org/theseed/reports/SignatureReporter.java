/**
 *
 */
package org.theseed.reports;

import java.io.PrintWriter;
import java.util.Map;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * This is the base class for all reports produced by the SignatureProcessor.
 *
 * @author Bruce Parrello
 *
 */
public abstract class SignatureReporter {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(SignatureReporter.class);
    /** print writer for text reports */
    protected PrintWriter writer;
    /** controlling command processor */
    protected IParms processor;
    /** map of signature IDs to names */
    private Map<String, String> nameMap;
    /** number of genomes in IN group */
    private int total1;
    /** number of genomes in OUT group */
    private int total2;

    /**
     * Interface for controlling processors, used to get tuning parameters
     */
    public interface IParms {

    }

    /**
     * Enum for report types
     */
    public static enum Type {
        COUNTS {
            @Override
            public SignatureReporter create(PrintWriter writer, IParms processor) {
                return new CountSignatureReporter(writer, processor);
            }
        };

        /**
         * @return a signature reporting object
         *
         * @param writer		output print writer
         * @param processor		controlling command processor
         */
        public abstract SignatureReporter create(PrintWriter writer, IParms processor);
    }

    /**
     * Construct a signature reporting object.
     *
     * @param writer		print writer for output
     * @param processor		controlling command processor
     */
    public SignatureReporter(PrintWriter writer, IParms processor) {
        this.writer = writer;
        this.processor = processor;
    }

    /**
     * Initialize the report.
     *
     * @param sigMap	map of signature class IDs to names
     * @param size1		number of IN-group genomes
     * @param size2		number of OUT-group genomes
     */
    public void openReport(Map<String, String> sigMap, int size1, int size2) {
        this.nameMap = sigMap;
        this.total1 = size1;
        this.total2 = size2;
        this.initReport();
    }

    /**
     * Initialize the report.  This usually involves writing the header.
     */
    protected abstract void initReport();

    /**
     * Output a single signature class.
     *
     * @param signature		signature class to display
     * @param count1		number of occurrences in IN group
     * @param count2		number of occurrences in OUT group
     */
    public abstract void showClass(String signature, int count1, int count2);

    /**
     * Complete the report.
     */
    public abstract void closeReport();

    /**
     * Skip a line to prepare for a secondary report.
     */
    public void space() {
        this.writer.println();
    }

    /**
     * @return the name of a signature class
     *
     * @param signature		class of interest
     */
    public String getName(String signature) {
        return this.nameMap.getOrDefault(signature, "?");
    }

    /**
     * @return the specified count as a percentage of the IN group
     *
     * @param count		count to convert
     */
    public double getPercent1(int count) {
        return getPercent(count, this.total1);
    }

    /**
     * @return the specified count as a percentage of the OUT group
     *
     * @param count		count to convert
     */
    public double getPercent2(int count) {
        return getPercent(count, this.total2);
    }

    /**
     * @return the count as a safe percentage of the total
     *
     * @param count		numerator
     * @param total		denominator
     */
    private static double getPercent(int count, int total) {
        double retVal = 100.0;
        if (total > 0)
            retVal = count * 100.0 / total;
        return retVal;
    }

}
