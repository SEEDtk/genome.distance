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
 * The input file is tab-delimited with headers and consists of two columns-- a method type and a parameter
 * string.  The parameter string is free-form and dependent on the method type; generally it consists of
 * space-delimited values.
 *
 * The positional parameters are the name of the genome source file or directory and the name of a role
 * definition file for role-parsing.  The command-line options are as follows:
 *
 * -h	display command-line usage
 * -v	display more frequent progress messages
 * -i	input file name (if not STDIN)
 * -o	output file name (if not STDOUT)
 *
 *
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
