/**
 *
 */
package org.theseed.genome.distance.methods;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import org.theseed.counters.CountMap;

/**
 * This object manages a list of genome ID pairs, and tracks the number of times each genome is used.
 * The pairs are internally sorted by most common genome to least, and then the pairs themselves are
 * sorted by the same criterion on the first genome.  The idea is to create an ordering that decreases
 * the number of genome loads required.
 *
 * @author Bruce Parrello
 *
 */
public class GenomePairList implements Iterable<GenomePairList.Pair> {

    // FIELDS
    /** list of pair objects */
    private List<Pair> pairs;
    /** genome occurrence counts */
    private CountMap<String> counts;
    /** set of pairs already inserted */
    private Set<Pair> pairSet;
    /** default list size */
    private static final int LIST_SIZE = 1000;

    /**
     * Object representing a pair of genome IDs.
     */
    public class Pair implements Comparable<Pair> {

        /** first genome ID */
        private String id1;
        /** second genome ID */
        private String id2;

        /**
         * Construct a genome pair.  Pairs are unordered, so we initialize them
         * to lexically sorted.  This insures that they remain lexically sorted
         * unless the second one is strictly more common than the first.
         *
         * @param	gId1	first genome ID
         * @param	gId2	second genome ID
         */
        public Pair(String gId1, String gId2) {
            if (gId1.compareTo(gId2) > 0) {
                this.id1 = gId2;
                this.id2 = gId1;
            } else {
                this.id1 = gId1;
                this.id2 = gId2;
            }
            // Count the genomes.
            GenomePairList.this.counts.count(gId1);
            GenomePairList.this.counts.count(gId2);
        }

        /**
         * Fix the genome pair to put the most common genome first.
         */
        public void fixup() {
            if (GenomePairList.this.counts.getCount(this.id1) < GenomePairList.this.counts.getCount(this.id2)) {
                String buffer = this.id1;
                this.id1 = this.id2;
                this.id2 = buffer;
            }
        }

        /**
         * Compare two pairs.  We want the pair whose first genome is most common to be first, and after
         * that, we sort by the genome IDs themselves.
         */
        @Override
        public int compareTo(Pair o) {
            int count1 = GenomePairList.this.counts.getCount(this.id1);
            int count1o = GenomePairList.this.counts.getCount(o.id1);
            int retVal = (count1o - count1);
            if (retVal == 0) {
                retVal = this.id1.compareTo(o.id1);
                if (retVal == 0)
                    retVal = this.id2.compareTo(o.id2);
            }
            return retVal;
        }

        /**
         * @return the first (most common) genome ID
         */
        public String getId1() {
            return this.id1;
        }

        /**
         * @return the second (less common) genome ID
         */
        public String getId2() {
            return this.id2;
        }

        @Override
        public int hashCode() {
            final int prime = 31;
            int result = 1;
            result = prime * result + ((this.id1 == null) ? 0 : this.id1.hashCode());
            result = prime * result + ((this.id2 == null) ? 0 : this.id2.hashCode());
            return result;
        }

        @Override
        public boolean equals(Object obj) {
            if (this == obj) {
                return true;
            }
            if (!(obj instanceof Pair)) {
                return false;
            }
            Pair other = (Pair) obj;
            if (this.id1 == null) {
                if (other.id1 != null) {
                    return false;
                }
            } else if (!this.id1.equals(other.id1)) {
                return false;
            }
            if (this.id2 == null) {
                if (other.id2 != null) {
                    return false;
                }
            } else if (!this.id2.equals(other.id2)) {
                return false;
            }
            return true;
        }

    }

    /**
     * Create an empty genome pair list.
     */
    public GenomePairList() {
        this.counts = new CountMap<String>();
        this.pairSet = new HashSet<Pair>(LIST_SIZE);
        this.pairs = null;
    }

    /**
     * Add a pair to the pair list.
     *
     * @param id1	ID of first genome
     * @param id2	ID of second genome
     */
    public void addPair(String id1, String id2) {
        Pair newPair = this.new Pair(id1, id2);
        this.pairSet.add(newPair);
    }

    /**
     * Set up the pair list for iteration.
     */
    public void prepare() {
        // Transfer the pairs from the set to the list.  The set helped us to remove duplicates.
        this.pairs = new ArrayList<Pair>(this.pairSet);
        // Free up the pair-set memory.
        this.pairSet = null;
        // Fix up the pairs using the finished counts.
        this.pairs.stream().forEach(x -> x.fixup());
        // Sort the pair list.
        Collections.sort(this.pairs);
    }

    @Override
    public Iterator<Pair> iterator() {
        return this.pairs.iterator();
    }

    /**
     * @return the number of pairs in this list
     */
    public int size() {
        int retVal = 0;
        if (this.pairSet != null)
            retVal = this.pairSet.size();
        else
            retVal = this.pairs.size();
        return retVal;
    }

    /**
     * @return the genome pair at the specified list position
     *
     * @param i		position to query
     */
    public Pair get(int i) {
        return this.pairs.get(i);
    }

    /**
     * @return a set of the IDs in this pair list
     */
    public SortedSet<String> getIdSet() {
        SortedSet<String> retVal = new TreeSet<String>();
        // We have to figure out where the pairs are.
        Collection<Pair> pairings = this.pairSet;
        if (pairings == null) pairings = this.pairs;
        // Loop through the pairs.
        for (var pair : pairings) {
            retVal.add(pair.id1);
            retVal.add(pair.id2);
        }
        return retVal;
    }

}
