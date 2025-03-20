/**
 *
 */
package org.theseed.genome.distance;

import java.util.HashMap;
import java.util.Map;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.theseed.io.TabbedLineReader;

/**
 * This object describes a group type for the distance check processor. For each group type,
 * it contains the input column index, the group type name, a map of genome IDs to group IDs,
 * and the in-group and out-group distance statistics.
 *
 * @author Bruce Parrello
 *
 */
public class GroupTypeSpec {

	// FIELDS
	/** input file column index for group ID */
	private int colIdx;
	/** group type name */
	private String typeName;
	/** map of genome IDs to group IDs */
	private Map<String, String> genomeMap;
	/** in-group distance statistics */
	private SummaryStatistics inStats;
	/** in-group ones */
	private int inOnes;
	/** out-group distance statistics */
	private SummaryStatistics outStats;
	/** out-group ones */
	private int outOnes;
	/** number of genomes without a group */
	private int badPairCount;

	/**
	 * Create a new group-type specification with no data in it.
	 *
	 * @param col		input column index
	 * @param name		group type name
	 */
	public GroupTypeSpec(int col, String name) {
		this.colIdx = col;
		this.typeName = name;
		// Create the genome map and the statistics objects.
		this.genomeMap = new HashMap<String, String>();
		this.inStats = new SummaryStatistics();
		this.outStats = new SummaryStatistics();
		// Clear the counters.
		this.badPairCount = 0;
		this.inOnes = 0;
		this.outOnes = 0;
	}

	/**
	 * Process a genome input line. We use this to update the genome map.
	 *
	 * @param line		input file line with the genome ID in column 1
	 */
	public void addGenome(TabbedLineReader.Line line) {
		String genomeId = line.get(0);
		String groupId = line.get(this.colIdx);
		genomeMap.put(genomeId, groupId);
	}

	/**
	 * Process a distance between two genomes. This updates the statistics.
	 *
	 * @param genome1	ID of the first genome
	 * @param genome2	ID of the second genome
	 * @param dist		distance between the genomes
	 */
	public void recordDistance(String genome1, String genome2, double dist) {
		String group1 = this.genomeMap.get(genome1);
		String group2 = this.genomeMap.get(genome2);
		if (group1 == null || group2 == null)
			this.badPairCount++;
		else if (group1.equals(group2)) {
			// Here we are in-group.
			if (dist == 1.0)
				this.inOnes++;
			else
				this.inStats.addValue(dist);
		} else {
			// Here we are out-group.
			if (dist == 1.0)
				this.outOnes++;
			else
				this.outStats.addValue(dist);
		}
	}

	/**
	 * Erase the current statistics.
	 */
	public void clear() {
		this.inStats.clear();
		this.inOnes = 0;
		this.outStats.clear();
		this.outOnes = 0;
	}

	/**
	 * @return the number of in-group unit distances
	 */
	public int getInOnes() {
		return this.inOnes;
	}

	/**
	 * @return the number of out-group unit distances
	 */
	public int getOutOnes() {
		return this.outOnes;
	}

	/**
	 * @return the group type name
	 */
	public String getTypeName() {
		return this.typeName;
	}

	/**
	 * @return the in-group statistics
	 */
	public SummaryStatistics getInStats() {
		return this.inStats;
	}

	/**
	 * @return the out-group statistics
	 */
	public SummaryStatistics getOutStats() {
		return this.outStats;
	}

	/**
	 * @return the number of bad genome pairs
	 */
	public int getBadPairCount() {
		return this.badPairCount;
	}

}
