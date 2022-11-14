/**
 *
 */
package org.theseed.genome.distance;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Map;
import java.util.TreeMap;
import org.apache.commons.lang3.StringUtils;
import org.kohsuke.args4j.Argument;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.io.TabbedLineReader;
import org.theseed.utils.BasePipeProcessor;
import org.theseed.utils.ParseFailureException;
import org.theseed.utils.StringPair;

/**
 * This command merges two files with identical formats that have two key columns.  The "new" file is read on the standard
 * input.  The "old" file name is specified as a positional parameter.  Any input line in the old file with the same keys as
 * the new file is removed, and the remaining lines added to the end of the output file after all the new-file lines are
 * written, in lexical order.
 *
 * The positional parameters are the old file name and the names of the two key columns.  The keys are treated as an unordered
 * string pair, so A B will match B A.
 *
 * The command-line options are as follows;
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -i	input (new) file name (if not STDIN)
 * -o	output file name (if not STDOUT)
 *
 * @author Bruce Parrello
 *
 */
public class PairMergeProcessor extends BasePipeProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(PairMergeProcessor.class);
    /** array of labels from the old file */
    private String[] oldLabels;
    /** index of the first key column */
    private int key1idx;
    /** index of the second key column */
    private int key2idx;
    /** map of old file key pairs to output lines */
    private Map<StringPair, TabbedLineReader.Line> oldFileMap;

    // COMMAND-LINE OPTIONS

    /** name of the old input file */
    @Argument(index = 0, metaVar = "oldFile.tbl", usage = "name of the old input file (duplicate lines will be removed before merge)",
            required = true)
    private File oldFile;

    /** index (1-based) or name of the first key field */
    @Argument(index = 1, metaVar = "key1Col", usage = "index (1-based) or name of the first key field", required = true)
    private String key1Name;

    /** index (1-based) or name of the first key field */
    @Argument(index = 2, metaVar = "key2Col", usage = "index (1-based) or name of the second key field", required = true)
    private String key2Name;

    @Override
    protected void setPipeDefaults() {
    }

    @Override
    protected void validatePipeInput(TabbedLineReader inputStream) throws IOException {
        // Verify that the labels match.
        if (! Arrays.equals(inputStream.getLabels(), this.oldLabels))
            throw new IOException("Input stream columns do not match old input file " + this.oldFile + ".");
    }

    @Override
    protected void validatePipeParms() throws IOException, ParseFailureException {
        // First, make sure we can read the old input file.
        if (! this.oldFile.canRead())
            throw new FileNotFoundException("Old input file " + this.oldFile + " is not found or unreadable.");
        // We process the old file in its entirety here.
        try (TabbedLineReader oldStream = new TabbedLineReader(this.oldFile)) {
            log.info("Verifying old input file {}.", this.oldFile);
            // Verify the key columns.
            this.key1idx = oldStream.findField(this.key1Name);
            this.key2idx = oldStream.findField(this.key2Name);
            // Save the labels.
            this.oldLabels = oldStream.getLabels();
            // Now read in the whole file and store the string pairs.
            this.oldFileMap = new TreeMap<StringPair, TabbedLineReader.Line>();
            oldStream.stream().forEach(x -> this.oldFileMap.put(this.getLineKey(x), x));
        }
        log.info("{} unique key pairs found in old input file.", this.oldFileMap.size());
    }

    @Override
    protected void runPipeline(TabbedLineReader inputStream, PrintWriter writer) throws Exception {
        // Write the headers.
        writer.println(StringUtils.join(this.oldLabels, '\t'));
        // Get some counters.
        int lineCount = 0;
        int deleteCount = 0;
        // Loop through the new input file.
        log.info("Scanning new input records.");
        for (var line : inputStream) {
            StringPair newKey = this.getLineKey(line);
            if (this.oldFileMap.containsKey(newKey)) {
                this.oldFileMap.remove(newKey);
                deleteCount++;
            }
            // Write the line from the new file.
            writer.println(line.toString());
            lineCount++;
            if (log.isInfoEnabled() && lineCount % 5000 == 0)
                log.info("{} input lines copied.", lineCount);
        }
        // Now flush out the old-file lines.
        log.info("{} old-file lines deleted.  {} remaining.", deleteCount, this.oldFileMap.size());
        for (var line : this.oldFileMap.values())
            writer.println(line.toString());
    }

    /**
     * Extract the key fields from an input line in either file.
     *
     * @param line		input line of interest
     *
     * @return a StringPair of the key fields
     */
    private StringPair getLineKey(TabbedLineReader.Line line) {
        StringPair retVal = new StringPair(line.get(this.key1idx), line.get(this.key2idx));
        return retVal;
    }

}
