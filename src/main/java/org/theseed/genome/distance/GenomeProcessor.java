/**
 *
 */
package org.theseed.genome.distance;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.IntStream;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.BaseReportProcessor;
import org.theseed.basic.ParseFailureException;
import org.theseed.genome.Genome;
import org.theseed.genome.iterator.GenomeSource;
import org.theseed.sequence.GenomeKmers;

/**
 * This class uses kmer distance to compare all the genomes in one or more genome sources using dna
 * contig kmers.  This is much more memory-intensive than the protein-based comparison. The positional
 * parameters are the name of the directory containing the base genomes and the name of the genome directory
 * containing the comparison genomes.  Multiple comparison genome directories can be specified.
 *
 * The command-line options are as follows.
 *
 * -v	show more detail on the log
 * -h	display usage information
 * -K	kmer size to use
 * -m	maximum distance; if a comparison results in a distance greater than this value, it will not be output;
 * 		the default is 1.0, which outputs everything
 * -t	type of genome source (default DIR)
 *
 * @author Bruce Parrello
 *
 */
public class GenomeProcessor extends BaseReportProcessor {

    // FIELDS

    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(GenomeProcessor.class);
    /** list of kmers for main genomes */
    private List<GenomeKmers> mainKmers;

    // COMMAND-LINE OPTIONS

    /** kmer size */
    @Option(name="--kmerSize", aliases = { "-K", "--kmer" }, metaVar = "12", usage = "DNA kmer size")
    private int kmerSize;

    /** maximum distance to allow */
    @Option(name = "--maxDist", aliases = { "-m", "--max", "--distance" }, metaVar = "0.75",
            usage = "maximum acceptable distance for a neighboring genome")
    private double maxDist;

    /** genome source type */
    @Option(name = "--type", aliases = { "-t" }, usage = "genome source type")
    private GenomeSource.Type sourceType;

    /** name of the file or directory for the base genome source */
    @Argument(index=0, metaVar="gtoDir", required=true, usage="base genome source")
    private File baseDir;

    /** name of the files or directories for the other genome sources */
    @Argument(index=1, metaVar="gtoDir1 gtoDir2 ...", required=true, usage="directory of input GTOs")
    private List<File> genomeDirs;

    @Override
    protected void setReporterDefaults() {
        this.kmerSize = 21;
        this.maxDist = 0.9;
        this.sourceType = GenomeSource.Type.DIR;
    }

	@Override
	protected void validateReporterParms() throws IOException, ParseFailureException {
		// Verify the kmer size.
		if (this.kmerSize < 4)
			throw new ParseFailureException("Kmer size cannot be less than 4.");
		GenomeKmers.setKmerSize(this.kmerSize);
		log.info("Chosen kmer size is {}.", this.kmerSize);
		// Verify the maximum distance.
		if (this.maxDist <= 0.0 || this.maxDist > 1.0)
			throw new ParseFailureException("Maximum distance must be > 0 and <= 1.");
		// Validate the main genome source.
		if (! this.baseDir.exists())
			throw new FileNotFoundException("Main genome source \"" + this.baseDir + "\" is not found.");
		// Validate the alternate genome sources.
		for (File genomeDir : genomeDirs) {
			if (! genomeDir.exists())
				throw new FileNotFoundException("Genome source \"" + genomeDir + "\" is not found.");
		}
		// Now we begin the laborious process of loading the base genome kmers.
		try {
			GenomeSource baseGenomes = this.sourceType.create(this.baseDir);
			final int nGenomes = baseGenomes.size();
			int count = 0;
			this.mainKmers = new ArrayList<GenomeKmers>(nGenomes);
			log.info("Loading {} genomes from {}.", nGenomes, this.baseDir);
			for (Genome genome : baseGenomes) {
				count++;
				log.info("Processing genome {} of {}: {}.", count, nGenomes, genome);
				GenomeKmers kmers = new GenomeKmers(genome);
				this.mainKmers.add(kmers);
			}
		} catch (Exception e) {
			// Convert the missing-algorithm exceptions thrown by the MD5 stuff in GenomeKmers.
			throw new ParseFailureException(e.toString());
		}
	}

	@Override
	protected void runReporter(PrintWriter writer) throws Exception {
		// Write the output headers.
		writer.println("genome1\tgenome2\tdistance");
		int compares = 0;
		// Save the number main genomes and create the output array for distances.
		final int nMain = this.mainKmers.size();
		double[] distances = new double[nMain];
		// The basic strategy is to go through each secondary genome one at a time, and compare it to all
		// the base genomes in parallel.
		final int nDirs = this.genomeDirs.size();
		for (int idx = 0; idx < nDirs; idx++) {
			File dir = this.genomeDirs.get(idx);
			log.info("Loading genome directory {}.", dir);
			GenomeSource genomes = this.sourceType.create(dir);
			final int nGenomes = genomes.size();
			int gCount = 0;
			for (Genome genome : genomes) {
				gCount++;
				log.info("Processing {} genome {} of {}: {}.", dir, gCount, nGenomes, genome);
				// Compute all the distances for this genome.
				GenomeKmers kmers = new GenomeKmers(genome);
				IntStream.range(0, nMain).parallel().forEach(i -> distances[i] = kmers.distance(this.mainKmers.get(i)));
				// Output the results.
				String genome_id = genome.getId();
				for (int i = 0; i < nMain; i++) {
					writer.println(genome_id + "\t" + this.mainKmers.get(i).getGenomeId() + "\t" + distances[i]);
					compares++;
				}
			}
		}
		log.info("{} comparisons output.", compares);
	}

}
