/**
 *
 */
package org.theseed.genome.distance;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.lang3.StringUtils;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.BaseReportProcessor;
import org.theseed.basic.ParseFailureException;
import org.theseed.io.TabbedLineReader;

/**
 * This command processes output files from genome distance comparisons. Each output file is presumed to be for
 * a different type of comparison. These comparisons are then compared to in-group and out-group statistics
 * for various types of groups.
 *
 * There is a single tab-delimited input file (with headers) containing genome-grouping information. For each
 * grouping type (there can be more than one), there is a single column containing group IDs. The column
 * title will be the grouping name. Our first job will be to create maps from each genome to its groups of each
 * type. Then for each grouping type, we track in-group and out-group distances. The output report will report
 * the minimum, mean, maximum, and deviation edges of in-group and out-group for each grouping type within each
 * grouping file. The intent is to show which distance is best at predicting group membership.
 *
 * The positional parameters are the name of the grouping input file and the names of the distance files. Each
 * distance file must have genome IDs in the first two columns and distances in the third. (The columns to be
 * used for the grouping types are specified by the --cols command-line option. The default is just "3".) If
 * a distance file is a directory, all the "*.tbl" files in the directory will be parsed.
 *
 * The command-line options are as follows:
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -o	output file for report (if not STDOUT)
 *
 * --cols	a command-delimited list of columns, one per grouping type, containing the group IDs; each should
 * 			be the index (1-based) or name of the appropriate column
 *
 * @author Bruce Parrello
 *
 */
public class DistanceCheckProcessor extends BaseReportProcessor {

	// FIELDS
	/** logging facility */
	protected static Logger log = LoggerFactory.getLogger(DistanceCheckProcessor.class);
	/** array of grouping type specifications */
	private GroupTypeSpec[] groupSpecs;
	/** final set of distance file names */
	private List<File> distFiles;
	/** file filter for "tbl" files */
	private static final FilenameFilter TBL_FILTER = new FilenameFilter() {
		@Override
		public boolean accept(File dir, String name) {
			return name.endsWith(".tbl");
		}
	};


	// COMMAND-LINE OPTIONS

	/** comma-delimited list of column specifications */
	@Option(name = "--cols", metaVar = "genus,species",
			usage = "comma-delimited list of column specs (1-based index or name) for grouping columns")
	private String colSpecString;

	/** name of the grouping input file */
	@Argument(index = 0, metaVar = "genomes.tbl", usage = "name of input file containing genome IDs and groupings",
			required = true)
	private File genomeFile;

	/** names of the distance files */
	@Argument(index = 1, metaVar = "dists1.tbl dists2.tbl ...", usage = "names of the distance files/directories",
			required = true)
	private List<File> distFileSpecs;

	@Override
	protected void setReporterDefaults() {
		this.colSpecString = "3";
		this.distFileSpecs = new ArrayList<File>();
	}

	@Override
	protected void validateReporterParms() throws IOException, ParseFailureException {
		// Form a final list of the distance files.
		this.distFiles = new ArrayList<File>(this.distFileSpecs.size());
		for (File distFile : this.distFileSpecs) {
			if (distFile.isDirectory()) {
				File[] subFiles = distFile.listFiles(TBL_FILTER);
				for (File subFile : subFiles)
					this.checkInput(subFile);
			} else
				this.checkInput(distFile);
		}
		log.info("{} distance files found.", this.distFiles.size());
		// Now we need to load the genome IDs and groups. First we validate and open the genome input file.
		if (! this.genomeFile.canRead())
			throw new FileNotFoundException("Genome input file " + this.genomeFile + " is not found or unreadable.");
		try (TabbedLineReader gStream = new TabbedLineReader(this.genomeFile)) {
			// Break the column specifications into pieces.
			String[] colSpecs = StringUtils.split(colSpecString, ',');
			// Now that we know the number of grouping types, we can build the grouping-type array.
			this.groupSpecs = new GroupTypeSpec[colSpecs.length];
			// Loop through the column specifications, creating the group-type specifications.
			String[] labels = gStream.getLabels();
			for (int i = 0; i < colSpecs.length; i++) {
				String colSpec = colSpecs[i];
				int colIdx = gStream.findField(colSpec);
				this.groupSpecs[i] = new GroupTypeSpec(colIdx, labels[colIdx]);
			}
			// Now that we have all the group-type specifications built, we can process the
			// groupings.
			int gCount = 0;
			for (var line : gStream) {
				gCount++;
				for (var groupSpec : this.groupSpecs)
					groupSpec.addGenome(line);
			}
			log.info("Groupings stored for {} genomes.", gCount);
		}
	}

	/**
	 * Check an input distance file. If it's valid, add it to the input list.
	 *
	 * @param subFile	distance file to check
	 *
	 * @throws FileNotFoundException
	 */
	protected void checkInput(File subFile) throws FileNotFoundException {
		if (! subFile.canRead())
			throw new FileNotFoundException("Input distance file " + subFile + " is not found or unreadable.");
		this.distFiles.add(subFile);
	}

	@Override
	protected void runReporter(PrintWriter writer) throws Exception {
		// Now we process each file individually. For each grouping type, we will output in-group and out-group stats.
		writer.println("dist_file\tgroup_type\tin_out\tmin\tlow\tmean\thigh\tmax\tones");
		int fileCount = 0;
		for (File distFile : this.distFiles) {
			fileCount++;
			log.info("Processing file {} of {}: {}.", fileCount, this.distFiles.size(), distFile);
			long lastMsg = System.currentTimeMillis();
			// We use the file's base name for the report.
			String fileName = distFile.getName();
			try (TabbedLineReader distStream = new TabbedLineReader(distFile)) {
				// Clear the group specifications for the new file.
				Arrays.stream(this.groupSpecs).forEach(x -> x.clear());
				// Loop through the input file, recording distances.
				int lineCount = 0;
				for (var line : distStream) {
					lineCount++;
					String g1 = line.get(0);
					String g2 = line.get(1);
					double dist = line.getDouble(2);
					Arrays.stream(this.groupSpecs).forEach(x -> x.recordDistance(g1, g2, dist));
					long now = System.currentTimeMillis();
					if (now - lastMsg >= 5000) {
						log.info("{} distances read from {}.", lineCount, fileName);
						lastMsg = now;
					}
				}
				log.info("{} total distances read from {}.", lineCount, fileName);
				// Now we need to output the group data.
				for (GroupTypeSpec groupSpec : this.groupSpecs) {
					String groupType = groupSpec.getTypeName();
					this.writeStats(writer, fileName, groupType, "in", groupSpec.getInStats(), groupSpec.getInOnes());
					this.writeStats(writer, fileName, groupType, "out", groupSpec.getOutStats(), groupSpec.getOutOnes());
				}
			}
		}
		int allBad = Arrays.stream(this.groupSpecs).mapToInt(x -> x.getBadPairCount()).sum();
		log.info("{} bad pairs encountered.", allBad);
	}

	/**
	 * Write an output line for the report. We display the distance range for either an in-group or
	 * out-group of a given type in a given file.
	 *
	 * @param writer		output print wrtier
	 * @param fileName		label for the file of interest
	 * @param groupType		type of grouping
	 * @param inOut			"in" for the in-group or "out" for the out-group
	 * @param stats			statistics to display
	 * @param oneCount		number of unit distances
	 */
	private void writeStats(PrintWriter writer, String fileName, String groupType, String inOut,
			SummaryStatistics stats, int oneCount) {
		double min, max, mean, low, high;
		if (stats.getN() == 0) {
			// Here every single value was a unit distance.
			min = 1.0;
			max = 1.0;
			mean = 1.0;
			low = 1.0;
			high = 1.0;
		} else {
			// Get the min, max, and mean.
			min = stats.getMin();
			max = stats.getMax();
			mean = stats.getMean();
			// Use the standard deviation to compute the high and low probable extremes.
			double sDev = stats.getStandardDeviation();
			low = mean - sDev;
			high = mean + sDev;
		}
		// Write the results.
		writer.println(fileName + "\t" + groupType + "\t" + inOut + "\t"
				+ min + "\t" + low + "\t" + mean + "\t" + high + "\t" + max + "\t" + oneCount);
	}

}
