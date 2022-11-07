/**
 *
 */
package org.theseed.genome.distance;

import java.io.IOException;
import java.io.PrintWriter;

import org.theseed.io.TabbedLineReader;
import org.theseed.utils.BasePipeProcessor;
import org.theseed.utils.ParseFailureException;

/**
 * This command reads an input file describing genome-comparison methods and runs all the methods on pairs
 * of incoming genomes.  The output displays the results of each method on each genome pair.
 *
 * @author Bruce Parrello
 *
 */
public class MethodTableProcessor extends BasePipeProcessor {

    @Override
    protected void setPipeDefaults() {
        // TODO code for setPipeDefaults

    }

    @Override
    protected void validatePipeInput(TabbedLineReader inputStream) throws IOException {
        // TODO code for validatePipeInput

    }

    @Override
    protected void validatePipeParms() throws IOException, ParseFailureException {
        // TODO code for validatePipeParms

    }

    @Override
    protected void runPipeline(TabbedLineReader inputStream, PrintWriter writer) throws Exception {
        // TODO code for runPipeline

    }
    // FIELDS
    // TODO data members for MethodTableProcessor

    // TODO constructors and methods for MethodTableProcessor
}
