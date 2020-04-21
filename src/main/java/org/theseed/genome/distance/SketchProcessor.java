/**
 *
 */
package org.theseed.genome.distance;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.io.TabbedLineReader;
import org.theseed.sequence.SequenceKmers;
import org.theseed.sequence.hash.Bucket;
import org.theseed.sequence.hash.Sketch;

/**
 * This class will convert a file of proteins into sketches.  Each sketch is associated with a group ID, which could be
 * a role, a protein family ID, or a gene name.  The sketches will be output to a file that can be read back in.
 *
 * The positional parameter is the name of the output file.
 *
 * The input file should be tab-delimited, with group IDs in one column and protein sequences
 * in the other.
 *
 * The command-line options are as follows.
 *
 * -h	display command usage
 * -v	show more detailed progress messages
 * -K	protein kmer size; the default is 8
 * -i	input file containing protein families (default STDIN)
 * -c	index (1-based) or name of the input column containing group IDs
 * -p	index (1-based) or name of the input column containing protein sequences
 * -w 	sketch size (width); the default is 360
 *
 * @author Bruce Parrello
 *
 */
public class SketchProcessor extends ProteinKmerReader {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(SketchProcessor.class);

    // COMMAND-LINE OPTIONS

    /** sketch size */
    @Option(name = "-w", aliases = { "--width", "--sketchSize" }, metaVar = "400", usage = "sketch size for each protein")
    private int width;

    /** output file name */
    @Argument(index = 0, metaVar = "outFile.ser", usage = "output file name", required = true)
    private File outFile;


    @Override
    protected void setDefaults() {
        this.initProteinParms();
        this.width = 360;
    }

    @Override
    protected boolean validateParms() throws IOException {
        this.validateProteinParms();
        if (! this.outFile.exists()) {
            // Verify that we can create this file.
            this.outFile.createNewFile();
        } else if (! this.outFile.canWrite())
            throw new FileNotFoundException("Cannot write to output file " + this.outFile + ".");
        if (this.width < 10)
            throw new IllegalArgumentException("Sketch width cannot be less than 10.");
        return true;
    }

    @Override
    protected void processProteins() throws IOException {
        // Create the bucket.
        Bucket outBucket = new Bucket();
        // This will count the proteins read.
        int protCount = 0;
        // Fill the bucket with sketches.
        for (TabbedLineReader.Line line : this.input()) {
            SequenceKmers protKmers = this.getProteinKmers(line);
            String group = this.getGroupId(line);
            Sketch sketch = new Sketch(protKmers.hashSet(this.width), group);
            outBucket.add(sketch);
            protCount++;
            if (protCount % 500 == 0)
                log.info("{} proteins processed.", protCount);
        }
        // Write the sketches to the output.
        log.info("Writing {} sketches to {}.", protCount, this.outFile);
        outBucket.save(this.outFile);
        log.info("All done.");
    }


}
