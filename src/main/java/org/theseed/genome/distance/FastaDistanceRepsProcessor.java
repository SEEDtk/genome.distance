/**
 *
 */
package org.theseed.genome.distance;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.BaseReportProcessor;
import org.theseed.basic.ParseFailureException;
import org.theseed.sequence.FastaInputStream;
import org.theseed.sequence.KmerType;
import org.theseed.sequence.Sequence;
import org.theseed.sequence.SequenceKmers;

/**
 * This subcommand takes as input a simple FASTA file and outputs the sequences believed to
 * be representative based on kmer distance.
 *
 * The FASTA file should be presented on the standard input. The distances will be written
 * to the standard output.
 *
 * The command-line options are as follows:
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -o	output file for report (if not STDOUT)
 * -i	input FASTA file (if not STDIN)
 * -K	kmer size to use (default 21 for DNA/RNA, 8 for protein)
 * -d	maximum distance allowed from a representative (default 0.97)
 *
 * --type	sequence type (default DNA)
 *
 * @author Bruce Parrello
 *
 */
public class FastaDistanceRepsProcessor extends BaseReportProcessor {

	// FIELDS
	/** logging facility */
	protected static Logger log = LoggerFactory.getLogger(FastaDistanceProcessor.class);
	/** current map of representative sequence IDs to kmer hashes */
	private Map<String, SequenceKmers> repMap;

	// COMMAND-LINE OPTIONS

	/** input FASTA file (if not STDIN) */
	@Option(name = "--input", aliases = { "-i" }, usage = "input FASTA file (if not STDIN)")
	private File inFile;

	/** kmer size to use */
	@Option(name = "--kSize", aliases = {"--kmerSize", "-K" }, usage = "kmer size to use; 0 for sequence type default")
	private int kmerSize;

	/** maximum distance allowed from a representative to be considered a neighbor */
	@Option(name = "--dist", aliases = { "--maxDist", "-d" }, usage = "maximum distance a neighbor can be from a representative")
	private double maxDist;

	/** sequence type */
	@Option(name = "--type", usage = "input sequence type")
	private KmerType seqType;

	@Override
	protected void setReporterDefaults() {
		this.inFile = null;
		this.kmerSize = 0;
		this.seqType = KmerType.DNA;
		this.maxDist = 0.97;
	}

	@Override
	protected void validateReporterParms() throws IOException, ParseFailureException {
		// Process the kmer size default.
		if (this.kmerSize == 0)
			this.kmerSize = this.seqType.getKmerSize();
		// Validate the kmer size.
		if (this.kmerSize < 2)
			throw new ParseFailureException("Kmer size must be at least 2.");
		// Insure the input file (if any) is readable.
		if (this.inFile != null && ! this.inFile.canRead())
			throw new FileNotFoundException("Input file " + this.inFile + " is not found or invalid.");
		// Create repgen map.
		this.repMap = new HashMap<String, SequenceKmers>(100);
	}

	/**
	 * @return the appropriate FASTA input stream
	 *
	 * @throws IOException
	 */
	private FastaInputStream openInput() throws IOException {
		FastaInputStream retVal;
		if (this.inFile == null) {
			log.info("Reading {} sequences from the standard input.", this.seqType);
			retVal = new FastaInputStream(System.in);
		} else {
			log.info("Reading {} sequences from {}.", this.seqType, this.inFile);
			retVal = new FastaInputStream(this.inFile);
		}
		return retVal;
	}

	@Override
	protected void runReporter(PrintWriter writer) throws Exception {
		int pairCount = 0;
		int seqCount = 0;
		// Write the output header.
		writer.println("seq\tname");
		// Loop until we run out of sequences.
		try (FastaInputStream seqStream = this.openInput()) {
			long lastMsg = System.currentTimeMillis();
			for (Sequence seq : seqStream) {
				seqCount++;
				// Get the kmers for this sequence.
				SequenceKmers seqKmers = this.seqType.createKmers(seq.getSequence(), this.kmerSize);
				// Check for a representative.
				Iterator<SequenceKmers> iter = this.repMap.values().iterator();
				boolean repFound = false;
				while (! repFound && iter.hasNext()) {
					SequenceKmers repKmers = iter.next();
					double dist = repKmers.distance(seqKmers);
					if (dist <= this.maxDist) {
						// Here we've found a representative.
						repFound = true;
					}
					pairCount++;
					long now = System.currentTimeMillis();
					if (now - lastMsg >= 10000) {
						log.info("{} pairs checked, {} representatives found for {} sequences.",
								pairCount, this.repMap.size(), seqCount);
						lastMsg = now;
					}
				}
				if (! repFound) {
					// Here we have a new representative.
					writer.println(seq.getLabel() + "\t" + seq.getComment());
					repMap.put(seq.getLabel(), seqKmers);
				}
			}
		}
		log.info("{} representatives found for {} sequences.", this.repMap.size(), seqCount);
	}

}
