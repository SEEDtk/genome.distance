/**
 *
 */
package org.theseed.genome.distance;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Set;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Genome;
import org.theseed.genome.iterator.GenomeSource;
import org.theseed.io.TabbedLineReader;
import org.theseed.sequence.GenomeKmers;
import org.theseed.utils.BaseReportProcessor;
import org.theseed.utils.ParseFailureException;

/**
 * This command accepts as input a table of genome pairs and outputs the distances for each.  The
 * genomes themselves must be available from a genome source.  The genome pairs will be read from
 * the standard input, which should be a tab-delimited file with the genome IDs in the first two
 * columns.
 *
 * The positional parameter is the name of the genome source directory (or file).  The command-line
 * options are as follows.
 *
 * -h	display command-line usage
 * -v	show more frequent log messages
 * -o	name of the output file for the report (if not STDOUT)
 * -i	name of the input file containing the pairs (if not STDIN)
 * -t	type of genome source (PATRIC, MASTER, DIR); the default is DIR
 * -R	role definition file for protein comparisons
 * -K	kmer size for contig comparisons
 *
 * --skipMissing	if specified, missing genomes will be skipped; the default is to terminate the
 * 					process
 *
 * --method		comparison method; the default is PROTEIN
 *
 *
 * @author Bruce Parrello
 *
 */
public class DistanceTableProcessor extends BaseReportProcessor implements Measurer.IParms {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(DistanceTableProcessor.class);
    /** input file stream */
    private TabbedLineReader inStream;
    /** list of genome pairs found */
    private GenomePairList pairList;
    /** source for finding genomes */
    private GenomeSource genomes;
    /** measurer for current first-genome */
    private Measurer current;
    /** set of genome IDs in the genome source */
    private Set<String> genomesFound;

    // COMMAND-LINE OPTIONS

    /** input file (if not STDIN) */
    @Option(name = "--input", aliases = { "-i" }, metaVar = "pairFile.tbl", usage = "file of genome pairs to process")
    private File inFile;

    /** type of input genome source */
    @Option(name = "--type", aliases = { "-t" }, usage = "type of genome source")
    private GenomeSource.Type inType;

    /** input role definition file */
    @Option(name = "--roles", aliases = { "-R" }, usage = "role definition file (for method = PROTEIN")
    private File roleFile;

    /** DNA kmer size */
    @Option(name = "--dnaKmer", aliases = { "-K", "--kmer" }, usage = "DNA kmer size (for method = CONTIG")
    private int dnaKmerSize;

    /** method for computing distance */
    @Option(name = "--method", usage = "method for computing distance")
    private Measurer.Type method;

    /** TRUE if we want to skip missing genomes instead of throwing an error */
    @Option(name = "--skipMissing", usage = "if specified, missing genomes will be skipped instead of causing an error")
    private boolean skipFlag;

    /** source for genomes */
    @Argument(index = 0, metaVar = "inDir", usage = "source for genomes")
    private File inDir;

    @Override
    protected void setReporterDefaults() {
        this.dnaKmerSize = GenomeKmers.kmerSize();
        this.roleFile = null;
        this.inFile = null;
        this.method = Measurer.Type.PROTEIN;
        this.inType = GenomeSource.Type.DIR;
        this.skipFlag = false;
    }

    @Override
    protected void validateReporterParms() throws IOException, ParseFailureException {
        // Initialize the measurement engine.
        this.method.init(this);
        // Load the genome source.
        if (! this.inDir.exists())
            throw new FileNotFoundException("Genome source " + this.inDir + " does not exist.");
        this.genomes = this.inType.create(this.inDir);
        this.genomesFound = this.genomes.getIDs();
        // Verify the input stream.
        if (this.inFile == null) {
            log.info("Pairs will be read from the standard input.");
            this.inStream = new TabbedLineReader(System.in);
        } else {
            log.info("Pairs will be read from {}.", this.inFile);
            this.inStream = new TabbedLineReader(this.inFile);
        }
        // Create the genome pair list.
        this.pairList = new GenomePairList();
    }

    @Override
    protected void runReporter(PrintWriter writer) throws Exception {
        try {
            // First, we read all the pairs from the standard input.
            for (TabbedLineReader.Line line : this.inStream) {
                // Insure we have both genomes in this pair.
                String id1 = line.get(0);
                String id2 = line.get(1);
                boolean ok = (this.verify(id1) && this.verify(id2));
                // Queue the pair for later processing.
                if (ok)
                    this.pairList.addPair(id1, id2);
            }
            // Insure we have something to do.
            if (this.pairList.size() == 0)
                log.warn("No pairs found in input.");
            else {
                log.info("{} pairs found in input.", this.pairList.size());
                // Prepare the list for processing.
                this.pairList.prepare();
                // Now we want to run through the list.  We start by writing the output header.
                writer.println("id1\tname1\tid2\tname2\tsimilarity");
                // Cache the first genome for the first pair.  The pairs are sorted by first genome,
                // with the most common one first.
                GenomePairList.Pair pair0 = this.pairList.get(0);
                this.current = this.getMeasurer(pair0.getId1());
                // Loop through the pairs.
                for (GenomePairList.Pair pair : this.pairList) {
                    String id1 = pair.getId1();
                    // Insure we have the correct measurer in the cache.
                    if (! this.current.getId().contentEquals(id1))
                        this.current = this.getMeasurer(id1);
                    // Get the second genome.
                    Genome otherGenome = this.genomes.getGenome(pair.getId2());
                    // Now compare the genomes.
                    double similarity = this.current.computePercentSimilarity(otherGenome);
                    writer.format("%s\t%s\t%s\t%s\t%8.4f%n", this.current.getId(), this.current.getName(),
                            otherGenome.getId(), otherGenome.getName(), similarity);
                }
            }
        } finally {
            // Close the input stream.  We do this explicitly in case this command is part of a
            // larger command.
            this.inStream.close();
        }

    }

    /**
     * @return a measurement object for the specified genome
     *
     * @param id	ID of the genome of interest
     */
    private Measurer getMeasurer(String id) {
        Genome genome = this.genomes.getGenome(id);
        log.info("Loading measurer for {}.", genome);
        Measurer retVal = this.method.create(genome);
        return retVal;
    }

    /**
     * Verify that the specified genome is in the genome source.
     *
     * @param id	ID of the genome to verify
     *
     * @returns TRUE if the genome is OK, FALSE if it is missing.
     *
     * @throws IOException
     */
    private boolean verify(String id) throws IOException {
        boolean retVal = true;
        if (! this.genomesFound.contains(id)) {
            if (this.skipFlag) {
                retVal = false;
                log.warn("Genome {} is missing from source", id);
            } else
                throw new IOException(String.format("Genome %s not found in source %s.", id, this.inDir));
        }
        return retVal;
    }

    @Override
    public File getRoleFile() {
        return this.roleFile;
    }

    @Override
    public int getKmerSize() {
        return this.dnaKmerSize;
    }

}
