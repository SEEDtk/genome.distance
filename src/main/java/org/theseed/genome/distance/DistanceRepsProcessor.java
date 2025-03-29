/**
 *
 */
package org.theseed.genome.distance;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.time.Duration;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.ParseFailureException;
import org.theseed.counters.CountMap;
import org.theseed.genome.Genome;
import org.theseed.genome.iterator.GenomeSource;
import org.theseed.sequence.GenomeKmers;
import org.theseed.utils.BaseMultiReportProcessor;

/**
 * This command will be used to classify genomes into representative sets using genome distance. It will build
 * a list of genomes that are all distant from each other and then find which of those is closest to each input
 * genome.
 *
 * It is worth noting that this is intended for viruses, where genome distance is a more practical measure.
 *
 * The positional parameters are the names of the genome sources. Multiple reports will be output.
 *
 * The command-line options are as follows:
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -D	output directory name
 * -K	kmer size to use (default 9)
 * -t	type of genome source (default DIR)
 *
 * --clear	erase the output directory before processing
 * --dist	distance to use; a distance greater than this is considered too far for inclusion (default 0.97)
 *
 * @author Bruce Parrello
 *
 */
public class DistanceRepsProcessor extends BaseMultiReportProcessor {

	// FIELDS
	/** logging facility */
	protected static Logger log = LoggerFactory.getLogger(DistanceRepsProcessor.class);
	/** map of genome IDs to genome kmers for the representatives */
	private Map<String, GenomeKmers> repMap;
	/** list of genome sources */
	private List<GenomeSource> genomeSources;
	/** total number of genomes in all sources */
	private int gTotal;
	/** default result object */
	private final static Result NULL_RESULT = new Result();

	// COMMAND-LINE OPTIONS

	/** kmer size to use for distance measurement */
	@Option(name = "--kmerSize", aliases = { "-K", "--kmer" }, metaVar = "20",
			usage = "kmer size to use for distance computation")
	private int kmerSize;

	/** type of input genome sources */
	@Option(name = "--sourceType", aliases = { "--type", "-t" }, usage = "type of genome sources")
	private GenomeSource.Type sourceType;

	/** maximum distance to be included in a representative's neighborhood */
	@Option(name = "--dist", metaVar = "0.9", usage = "maximum distance for a representative neighborhood")
	private double maxDist;

	/** names of the genome sources */
	@Argument(index = 0, metaVar = "inDir1 inDir2 ...", usage = "file or directory names of the genome sources",
			required = true)
	private List<File> inDirs;

	/**
	 * This dinky class tracks a GenomeKmers object and a distance.
	 */
	protected static class Result {

		/** target representative */
		private GenomeKmers repGen;
		/** distance to source */
		private double distance;

		/**
		 * Compute the distance result for a genome.
		 *
		 * @param source	source genome's kmers
		 * @param repKmers	repgen genome's kmers
		 */
		protected Result(GenomeKmers source, GenomeKmers repKmers) {
			this.distance = source.distance(repKmers);
			this.repGen = repKmers;
		}

		/**
		 * Create a null result.
		 */
		protected Result() {
			this.distance = 1.0;
			this.repGen = null;
		}

		/**
		 * Merge two results and return the one with the minimum distance.
		 *
		 * @param other		other result to merge into this one
		 *
		 * @return the minimal-distance result
		 */
		protected Result merge(Result other) {
			return (this.distance <= other.distance ? this : other);
		}

		/**
		 * @return the repgen's kmer object
		 */
		public GenomeKmers getRepGen() {
			return this.repGen;
		}

		/**
		 * @return the distance
		 */
		public double getDistance() {
			return this.distance;
		}

	}

	@Override
	protected File setDefaultOutputDir(File curDir) {
		return new File(curDir, "repDb");
	}

	@Override
	protected void setMultiReportDefaults() {
		this.inDirs = new ArrayList<File>();
		this.maxDist = 0.97;
		this.kmerSize = 9;
		this.sourceType = GenomeSource.Type.DIR;
	}

	@Override
	protected void validateMultiReportParms() throws IOException, ParseFailureException {
		// Validate the kmer size and save it.
		if (this.kmerSize < 4)
			throw new ParseFailureException("Kmer size must be at least 4.");
		GenomeKmers.setKmerSize(this.kmerSize);
		// Validate the maximum distance.
		if (this.maxDist <= 0.0 || this.maxDist >= 1.0)
			throw new ParseFailureException("Distance must be strictly between 0 and 1.");
		// Now set up the genome sources.
		this.gTotal = 0;
		this.genomeSources = new ArrayList<GenomeSource>(this.inDirs.size());
		for (File inDir : this.inDirs) {
			if (! inDir.exists())
				throw new FileNotFoundException("Genome source " + inDir + " is not found.");
			GenomeSource genomes = this.sourceType.create(inDir);
			this.gTotal += genomes.size();
			log.info("{} genomes found in {}.", genomes.size(), inDir);
			this.genomeSources.add(genomes);
		}
		log.info("{} total genomes found in all sources.", this.gTotal);
		// Initialize the representative-genome map.
		this.repMap = new HashMap<String, GenomeKmers>(500);
	}

	@Override
	protected void runMultiReports() throws Exception {
		log.info("Starting first pass to find representatives.");
		int gCount = 0;
		long lastMsg = System.currentTimeMillis();
		long start = lastMsg;
		// During the first pass, we find the representatives.
		for (GenomeSource genomes : this.genomeSources) {
			for (Genome genome : genomes) {
				gCount++;
				GenomeKmers kmers = new GenomeKmers(genome);
				// Find out if genome belongs with a current representative.
				boolean belongs = this.repMap.values().parallelStream().anyMatch(x -> x.distance(kmers) <= this.maxDist);
				if (! belongs) {
					// It doesn't, so add it to the representative set.
					this.repMap.put(genome.getId(), kmers);
				}
				long now = System.currentTimeMillis();
				if (now - lastMsg >= 5000) {
					log.info("{} of {} genomes processed, {} representatives chosen.", gCount, this.gTotal, this.repMap.size());
					lastMsg = now;
				}
			}
		}
		if (log.isInfoEnabled()) {
			long now = System.currentTimeMillis();
			Duration interval = Duration.ofMillis(now - start);
			log.info("{} total representatives found for {} genomes in {}.", this.repMap.size(), this.gTotal, interval);
		}
		// Now we make the second pass. Here we only load the genomes that are not representatives, and we find out
		// where they belong. We will produce a list file showing where each genome belongs and a stats file showing
		// the size of each representative's neighborhood. This map will track the neighborhood size.
		CountMap<String> neighborCounts = new CountMap<String>();
		// Compute the names we want to give to the output file.
		String namePrefix = String.format("rep%.4f_K%d", this.maxDist, this.kmerSize);
		// We will write the assignments to the following output file.
		File listFile = this.getOutFile(namePrefix + ".list.tbl");
		try (PrintWriter writer = new PrintWriter(listFile)) {
			writer.println("genome_id\tgenome_name\trep_id\trep_name\tdistance");
			gCount = 0;
			log.info("Assigning {} genomes to {} representatives.", this.gTotal, this.repMap.size());
			lastMsg = System.currentTimeMillis();
			for (GenomeSource genomeSource : this.genomeSources) {
				Collection<String> genomeIDs = genomeSource.getIDs();
				for (String genomeID : genomeIDs) {
					gCount++;
					String genomeName;
					double dist;
					// Find out if this genome is a repgen.
					GenomeKmers repFound = this.repMap.get(genomeID);
					if (repFound != null) {
						// Yes, so it belongs in its own neighborhood. We can copy the name from the repgen map
						// without loading the genome.
						genomeName = repFound.getGenomeName();
						dist = 0.0;
					} else {
						// Load this genome. This is an ugly construct, but it keeps the genome from being stuck
						// in memory.
						GenomeKmers newKmers = new GenomeKmers(genomeSource.getGenome(genomeID));
						// Here we have to find the best representative.
						Result best = this.repMap.values().parallelStream().map(x -> new Result(newKmers, x))
								.reduce(NULL_RESULT, (x,y) -> x.merge(y));
						// We know because of the first pass that this will not be null.
						repFound = best.getRepGen();
						dist = best.getDistance();
						// Save the genome name.
						genomeName = newKmers.getGenomeName();
					}
					// Now we have the name and ID of the neighboring genome, the GenomeKmers object for the
					// representative (which contains its name and ID), and the distance between the neighbor
					// and its representative.
					final String repID = repFound.getGenomeId();
					writer.println(genomeID + "\t" + genomeName + "\t" + repID + "\t" + repFound.getGenomeName()
							+ "\t" + dist);
					// Update the counters.
					neighborCounts.count(repID);
					long now = System.currentTimeMillis();
					if (now - lastMsg >= 5000) {
						log.info("{} genomes placed.", gCount);
						lastMsg = now;
					}
				}
			}
			log.info("{} total genomes placed.", gCount);
		}
		// Now we run through the count map and output the statistics.
		File statFile = this.getOutFile(namePrefix + ".stats.tbl");
		try (PrintWriter writer = new PrintWriter(statFile)) {
			log.info("Writing statistics file {}.", statFile);
			writer.println("rep_id\trep_name\tsize");
			var sortedCounts = neighborCounts.sortedCounts();
			for (var count : sortedCounts) {
				String repID = count.getKey();
				String repName = this.repMap.get(repID).getGenomeName();
				writer.println(repID + "\t" + repName + "\t" + count.getCount());
			}
		}
	}

}
