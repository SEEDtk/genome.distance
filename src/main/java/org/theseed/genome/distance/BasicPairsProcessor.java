/**
 *
 */
package org.theseed.genome.distance;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.Set;
import java.util.TreeSet;

import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.ParseFailureException;
import org.theseed.io.TabbedLineReader;
import org.theseed.utils.BasePipeProcessor;

/**
 * This is a simple program that reads a list of genome IDs and produces a list of pairs suitable for
 * the MethodTableProcessor.  Every genome will be paired with every other genome exactly once.  To
 * get a list of all pairs of repgen genomes, use the appropriate stats file.
 *
 * The input list of genome IDs should be tab-delimited, with headers, and the genome IDs in the first column.
 * The output list will have the genome IDs in the first two columns, with column names of "genome1" and "genome2".
 * In the default case, the input is the standard input and the output is the standard output.
 *
 * The command-line options are as follows:
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -i	input file name (if not STDIN)
 * -o	output file name (if not STDOUT)
 * -c	index (1-based) or name of input genome ID column (default "1")
 *
 * @author Bruce Parrello
 *
 */
public class BasicPairsProcessor extends BasePipeProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(BasicPairsProcessor.class);
    /** list of genome IDs */
    private Set<String> idSet;
    /** index of input genome ID column */
    private int idColIdx;

    // COMMAND-LINE OPTIONS

    /** index (1-based) or name of genome ID input column */
    @Option(name = "--col", aliases = { "-c", "--column" }, metaVar = "genome_id", usage = "index (1-based) or name of genome ID input column")
    private String idCol;

    @Override
    protected void setPipeDefaults() {
        this.idCol = "1";
    }

    @Override
    protected void validatePipeInput(TabbedLineReader inputStream) throws IOException {
        this.idColIdx = inputStream.findField(this.idCol);
        if (! inputStream.hasNext())
            throw new IOException("No genome IDs found in input.");
    }

    @Override
    protected void validatePipeParms() throws IOException, ParseFailureException {
    }

    @Override
    protected void runPipeline(TabbedLineReader inputStream, PrintWriter writer) throws Exception {
        // Write the output headers.
        writer.println("genome1\tgenome2");
        // Initialize our counters.
        int linesIn = 0;
        int linesOut = 0;
        // We use a sorted set so that the output is predictable.
        this.idSet = new TreeSet<String>();
        // Loop through the input.
        for (var line : inputStream) {
            // Get the genome ID.
            linesIn++;
            String genomeID = line.get(this.idColIdx);
            // Pair it with all the previous genomes.
            for (String otherID : this.idSet) {
                writer.println(genomeID + "\t" + otherID);
                linesOut++;
            }
            // Add it to the set for future pairings.
            this.idSet.add(genomeID);
        }
        // All done.
        log.info("{} lines read, {} unique genome IDs, {} pairs output.", linesIn, this.idSet.size(), linesOut);
    }

}
