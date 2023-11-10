/**
 *
 */
package org.theseed.genome.distance.methods;

import java.io.File;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.StringUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.ParseFailureException;
import org.theseed.genome.Contig;
import org.theseed.genome.Genome;
import org.theseed.p3api.P3Genome;
import org.theseed.sequence.DnaInputStream;
import org.theseed.sequence.FastaOutputStream;
import org.theseed.sequence.Sequence;
import org.theseed.sequence.blast.BlastDB;
import org.theseed.sequence.blast.BlastHit;
import org.theseed.sequence.blast.BlastParms;
import org.theseed.sequence.blast.DnaBlastDB;

/**
 * This distance method uses orthologous average nucleotide identity to compute distance.
 *
 * OrthoANI chops up each genome's contigs into equal-length chunks and isolates the bidirectional best
 * blast hits.  The standard ANI score is then the sum of the identity ratios divided by the number of blast-hit pairs.
 * The distance is 1 - score.  We require the match to have a minimum query identity of 35%, though this is tunable.  No
 * truncated chunks are allowed, only full-length ones.
 *
 * The parameter keywords are as follows:
 *
 * chunk	chunk size to use (default 1020)
 * minI		minimum percent query identity (default 30)
 * minM		minimum percent match length (default 35)
 *
 * @author Bruce Parrello
 *
 */
public class AniDistanceMethod extends DistanceMethod {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(AniDistanceMethod.class);
    /** temporary directory for working files */
    private File tempDir;
    /** minimum identity */
    private double minIdentity;
    /** minimum match length */
    private int minMatchLen;
    /** chunk size */
    private int chunkSize;
    /** blast parameters */
    private BlastParms blastParms;

    /**
     * The analysis method breaks the genome's contigs into chunks and stores them in a FASTA file in
     * the temporary directory.  The first time a genome is compared, it will be automatically converted
     * to a BLAST database if necessary.
     */
    public class Analysis extends Measurer {

        /** name of FASTA file */
        private File chunkFasta;
        /** blast database of the FASTA file */
        private DnaBlastDB blastDb;
        /** genetic code of this genome */
        private int gc;
        /** number of chunks */
        private int size;

        /**
         * Analyze a genome for ANI measurements.
         *
         * @param genome	genome to parse
         */
        public Analysis(Genome genome) {
            super(genome);
            // Denote that we have no blast database.
            this.blastDb = null;
            // Create the output FASTA file.
            this.chunkFasta = new File(tempDir, genome.getId() + ".fna");
            // Save the genetic code.
            this.gc = genome.getGeneticCode();
            // Initialize the chunk counter.
            this.size = 0;
            // Create a chunk ID prefix.
            String prefix = "G" + genome.getId() + "_";
            // Now we write the chunks.
            try (FastaOutputStream outStream = new FastaOutputStream(this.chunkFasta)) {
                for (Contig contig : genome.getContigs()) {
                    final String sequence = contig.getSequence();
                    // Insure we only process whole chunks.
                    final int n = sequence.length() - AniDistanceMethod.this.chunkSize;
                    // Loop through the sequence.
                    for (int i = 0; i < n; i += AniDistanceMethod.this.chunkSize)
                        this.writeChunk(outStream, prefix, sequence, i);
                }
            } catch (IOException e) {
                throw new UncheckedIOException(e);
            }
        }

        /**
         * Write the chunk at the specified position to the output stream.
         *
         * @param outStream		FASTA output stream
         * @param prefix		prefix to use for the chunk label
         * @param dna			contig DNA sequence
         * @param i				starting position of chunk
         *
         * @throws IOException
         */
        private void writeChunk(FastaOutputStream outStream, String prefix, String dna, int i) throws IOException {
            String chunk = StringUtils.substring(dna, i, i + AniDistanceMethod.this.chunkSize);
            this.size++;
            Sequence seq = new Sequence(prefix + Integer.toString(this.size), "", chunk);
            outStream.write(seq);
        }

        /**
         * @return the blast database for this genome's chunk FASTA
         */
        public BlastDB getBlastDB() {
            if (this.blastDb == null) try {
                this.blastDb = DnaBlastDB.create(this.chunkFasta, gc);
            } catch (IOException e) {
                throw new UncheckedIOException(e);
            } catch (InterruptedException e) {
                throw new RuntimeException(e);
            }
            return this.blastDb;
        }

        /**
         * @return a DNA input stream for this genome's chunk FASTA
         */
        public DnaInputStream getDnaStream() {
            DnaInputStream retVal;
            try {
                retVal = new DnaInputStream(this.chunkFasta, this.gc);
            } catch (IOException e) {
                throw new UncheckedIOException(e);
            }
            return retVal;
        }

    }

    @Override
    protected Measurer setupGenome(Genome genome) {
        return this.new Analysis(genome);
    }

    @Override
    protected void parseParms(Map<String, String> keywords) throws ParseFailureException, IOException {
        // Get the tuning parameters.
        this.chunkSize = this.getIntValue(keywords, "chunk", 1020);
        this.minIdentity = this.getDoubleValue(keywords, "minI", 30.0);
        this.minMatchLen = (int) (this.getDoubleValue(keywords, "minM", 35.0) * this.chunkSize / 100);
        // Create the temporary directory.
        this.tempDir = Files.createTempDirectory("aniB").toFile();
        log.info("Temporary ANI directory is {}.", this.tempDir);
        // Create the blast parameters.
        this.blastParms = new BlastParms().task("blastn").dustOff().xdrop_gap(150.0).penalty(-1).reward(1).maxE(1e-15)
                .minPercent(this.minIdentity).maxPerQuery(1).minLen(this.minMatchLen);
    }

    @Override
    public double getDistance(Measurer measurer, Measurer other) {
        AniDistanceMethod.Analysis m1 = (AniDistanceMethod.Analysis) measurer;
        AniDistanceMethod.Analysis m2 = (AniDistanceMethod.Analysis) other;
        // In general, the first genome is the one that will be reused, so we make it the blast data base
        // and make the other one a file stream.  We insure the file stream is closed.
        log.info("Accessing query and subject.");
        BlastDB db = m1.getBlastDB();
        List<BlastHit> hits;
        try (DnaInputStream dna = m2.getDnaStream()) {
            // Get the blast hits.
            log.info("Invoking BLASTN.");
            hits = db.blast(dna, this.blastParms);
        }
        // We need to track the best hits for each query sequence and for each subject sequence.  When the
        // two match, we have a BBH.
        var qMap = new HashMap<String, BlastHit>(Math.min(hits.size(), m2.size) * 4 / 3);
        var sMap = new HashMap<String, BlastHit>(Math.min(hits.size(), m1.size) * 4 / 3);
        // Note that we will only see hits that match our minimum criteria, because it's all in the blast parms.
        for (var hit : hits) {
            String query = hit.getQueryId();
            String subject = hit.getSubjectId();
            // Insure each sequence has its best hit mapped.  Note that when we are storing into the subject
            // map, both hits will have the same subject, and when we are storing into the query map, both hits
            // will have the same query.
            qMap.compute(query, (k, v) -> (v == null ? hit : this.merge(hit, v)));
            sMap.compute(subject, (k, v) -> (v == null ? hit : this.merge(hit, v)));
        }
        // We need the sum of the identity fractions and the number of BBHs.
        int sum = 0;
        int base = 0;
        int mCount = 0;
        // Loop through the best query hits.
        for (var hit : qMap.values()) {
            // Is this a BBH?
            String query = hit.getQueryId();
            BlastHit reciprocal = sMap.get(hit.getSubjectId());
            if (reciprocal.getQueryId().contentEquals(query)) {
                // Yes.  Get the average identity fraction and count the hit.
                sum += (hit.getNumIdentical() + reciprocal.getNumIdentical());
                base += hit.getAlignLen() + reciprocal.getAlignLen();
                mCount++;
            }
        }
        log.info("{} bidirectional best hits found.", mCount);
        // Compute the distance.
        double retVal;
        if (sum == 0.0)
            retVal = 1.0;
        else
            retVal = 1.0 - ((double) sum) / base;
        return retVal;
    }

    /**
     * Compare two blast hits and return the one with the highest identity count.
     * If both are the same, we keep the second one (which is the oldest).
     *
     * @param hit1	first blast hit
     * @param hit2	second blast hit
     *
     * @return the best hit
     */
    private BlastHit merge(BlastHit hit1, BlastHit hit2) {
        BlastHit retVal;
        if (hit1.getNumIdentical() > hit2.getNumIdentical())
            retVal = hit1;
        else
            retVal = hit2;
        return retVal;
    }

    @Override
    public P3Genome.Details getDetailLevel() {
        return P3Genome.Details.CONTIGS;
    }

    @Override
    public String getName() {
        String retVal = String.format("ANI_chunk%d.I%02d", this.chunkSize, (int) (this.minIdentity));
        if (this.minMatchLen > 0) {
            int pctMatch = (this.minMatchLen * 100 / this.chunkSize);
            retVal += String.format(",M>%02d", (int) pctMatch);
        }
        return retVal;
    }

    @Override
    public void close() throws Exception {
        // Insure we erase the temporary directory.
        if (this.tempDir != null && this.tempDir.isDirectory())
            FileUtils.forceDelete(this.tempDir);
    }

}
