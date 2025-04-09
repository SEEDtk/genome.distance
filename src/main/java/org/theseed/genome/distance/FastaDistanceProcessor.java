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
 * This subcommand takes as input a simple FASTA file and outputs the kmer distance between
 * each pair of sequences. This can be expensive if the sequences are long and there are
 * a lot of them.
 *
 * We will keep all of the sequences in memory, but we will only cache a fixed number of kmer
 * hashes in memory at a time. The sequences are stored in a list. After all the cached
 * sequences have been compared to all subsequence sequences, they will be removed from the
 * list, so the cache always refers to the first N list entries, where N is the batch size
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
 * -b	number of sequence kmer structures to cache in memory; also the maxmimum number of
 * 		parallel processes (default 20)
 *
 * --type	sequence type (default DNA)
 *
 * @author Bruce Parrello
 *
 */
public class FastaDistanceProcessor extends BaseReportProcessor {

	// FIELDS
	/** logging facility */
	protected static Logger log = LoggerFactory.getLogger(FastaDistanceProcessor.class);
	/** list of input sequences */
	private List<Sequence> sequences;
	/** current kmer cache */
	private SequenceKmers[] kmerCache;
	/** number of pairs output */
	private int pairCount;
	/** output print writer for report */
	private PrintWriter reportWriter;

	// COMMAND-LINE OPTIONS

	/** input FASTA file (if not STDIN) */
	@Option(name = "--input", aliases = { "-i" }, usage = "input FASTA file (if not STDIN)")
	private File inFile;

	/** kmer size to use */
	@Option(name = "--kSize", aliases = {"--kmerSize", "-K" }, usage = "kmer size to use; 0 for sequence type default")
	private int kmerSize;

	/** sequence batch size */
	@Option(name = "--batch", aliases = { "-b" }, usage = "batch size for kmer cache and parallelism")
	private int batchSize;

	/** sequence type */
	@Option(name = "--type", usage = "input sequence type")
	private KmerType seqType;

	@Override
	protected void setReporterDefaults() {
		this.inFile = null;
		this.kmerSize = 0;
		this.batchSize = 20;
		this.seqType = KmerType.DNA;
	}

	@Override
	protected void validateReporterParms() throws IOException, ParseFailureException {
		// Process the kmer size default.
		if (this.kmerSize == 0)
			this.kmerSize = this.seqType.getKmerSize();
		// Validate the kmer size.
		if (this.kmerSize < 2)
			throw new ParseFailureException("Kmer size must be at least 2.");
		// Validate the batch size.
		if (this.batchSize < 1)
			throw new ParseFailureException("Batch size must be at least 1.");
		// Now we need to read in the sequences.
		try (FastaInputStream inStream = this.openInput()) {
			this.sequences = new ArrayList<Sequence>(this.batchSize * 4);
			for (Sequence seq : inStream)
				this.sequences.add(seq);
			log.info("{} sequences read from input.", this.sequences.size());
		}
		// Create the sequence-kmer cache.
		this.kmerCache = new SequenceKmers[this.batchSize];
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
		} else if (! this.inFile.canRead())
			throw new FileNotFoundException("Input file " + this.inFile + " is not found or unreadable.");
		else {
			log.info("Reading {} sequences from {}.", this.seqType, this.inFile);
			retVal = new FastaInputStream(this.inFile);
		}
		return retVal;
	}

	@Override
	protected void runReporter(PrintWriter writer) throws Exception {
		int batchCount = 0;
		this.pairCount = 0;
		this.reportWriter = writer;
		// Write the output header.
		writer.println("seq1\tname1\tseq2\tname2\tdistance");
		// Loop until we run out of sequences.
		while (! this.sequences.isEmpty()) {
			batchCount++;
			log.info("Processing batch {} with {} sequences remaining.", batchCount, this.sequences.size());
			// Compute the size of the batch. It might be short if this is the last batch.
			final int bSize;
			if (this.batchSize < this.sequences.size())
				bSize = this.batchSize;
			else
				bSize = this.sequences.size();
			// Cache all the sequences in the batch.
			for (int i = 0; i < bSize; i++) {
				Sequence seq = this.sequences.get(i);
				SequenceKmers kmers = this.seqType.createKmers(seq.getSequence(), this.kmerSize);
				this.kmerCache[i] = kmers;
			}
			log.info("{} sequences cached. Computing distances.", bSize);
			IntStream.range(0, bSize).parallel()
					.forEach(i -> this.computePairs(i, this.sequences.get(i), this.kmerCache[i]));
			log.info("{} pairs computed in {} batches.", this.pairCount, batchCount);
			// Delete the sequences from the list.
			this.sequences.subList(0, bSize).clear();
		}
		// Insure we know the report writer is gone.
		this.reportWriter = null;
	}

	/**
	 * Compute all distances to the specified sequence and write them to the output.
	 *
	 * @param idx		list index of input sequence
	 * @param seq		input sequence
	 * @param seqKmers	cached kmer object for the input sequence
	 */
	private void computePairs(int idx, Sequence seq, SequenceKmers seqKmers) {
		// We need to compute the distance for each sequence that follows this one.
		final int nSeqs = this.sequences.size();
		for (int jdx = idx + 1; jdx < nSeqs; jdx++) {
			// Get the comparator sequence.
			Sequence seq2 = this.sequences.get(jdx);
			SequenceKmers seqKmers2;
			if (jdx < this.kmerCache.length)
				seqKmers2 = this.kmerCache[jdx];
			else
				seqKmers2 = this.seqType.createKmers(seq2.getSequence(), this.kmerSize);
			// Compute the distance.
			double distance = seqKmers.distance(seqKmers2);
			// Write the report line. This requires synchronizing.
			synchronized (this) {
				this.reportWriter.println(seq.getLabel() + "\t" + seq.getComment() + "\t" + seq2.getLabel()
				+ "\t" + seq2.getComment() + "\t" + distance);
				this.pairCount++;
			}
		}
	}

}
