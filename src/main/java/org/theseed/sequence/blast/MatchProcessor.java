/**
 *
 */
package org.theseed.sequence.blast;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import org.apache.commons.io.FileUtils;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Contig;
import org.theseed.genome.Genome;
import org.theseed.locations.Location;
import org.theseed.proteins.DnaTranslator;
import org.theseed.reports.NaturalSort;
import org.theseed.sequence.DnaDataStream;
import org.theseed.sequence.DnaInputStream;
import org.theseed.sequence.FastaOutputStream;
import org.theseed.sequence.Sequence;
import org.theseed.sequence.SequenceDataStream;
import org.theseed.sequence.blast.BlastDB;
import org.theseed.sequence.blast.BlastHit;
import org.theseed.sequence.blast.BlastParms;
import org.theseed.sequence.blast.DnaBlastDB;
import org.theseed.utils.BaseProcessor;

/**
 * This command takes DNA from the input and BLASTS it against the contigs of a genome, then parses proteins
 * out of the DNA that has good hits.
 *
 * The positional parameter is the name of the genome file.  The DNA FASTA file should come in via the standard
 * input.
 *
 * The command-line options are as follows.
 *
 * -v	show more detail on the log
 * -h	display usage information
 * -i	name of the input file for the DNA (the default is STDIN)
 * -b	number of input sequences to process at a time (the default is 10)
 * -x	distance to extend the genome hit on either side (the default is 50)
 *
 * --maxE		maximum permissible e-value (the default is 1e-10)
 * --minPct		minimum percent coverage of an incoming DNA sequence for a match (the default is 95)
 * --tempDir	temporary directory for BLAST databases; the default is "Temp" in the current directory
 * --maxGap		maximum gap between adjacent hits when they are to be joined.  The default is 500.
 * --filter		file name of a protein BLAST database to use for filtering out the input sequences. The
 * 				default is not to do any filtering
 * --fasta		optional FASTA file in which to store output proteins
 * --upstream	upstream region to include with output proteins
 *
 * @author Bruce Parrello
 *
 */
public class MatchProcessor extends BaseProcessor {

    // FIELDS
    /** logging facility */
    protected static final Logger log = LoggerFactory.getLogger(MatchProcessor.class);
    /** DNA input stream */
    private DnaInputStream inStream;
    /** input genome */
    private Genome genome;
    /** filtering BLAST datanase */
    private BlastDB filterDb;
    /** BLAST parameters for the filtering query */
    private BlastParms filterParms;
    /** BLAST parameters for the main query */
    private BlastParms mainParms;
    /** output stream for FASTA of proteins found */
    private FastaOutputStream fastaOut;

    // COMMAND-LINE OPTION

    /** input file */
    @Option(name = "-i", aliases = { "--input" }, metaVar = "rna.fasta", usage = "input file (if not STDIN)")
    private File inFile;

    /** number of DNA sequences to submit in each BLAST call */
    @Option(name = "-b", aliases = { "--batchSize", "--batch" }, metaVar = "1",
            usage = "number of input sequences to submit to each BLAST call")
    private int batchSize;

    /** distance to extend the genome hit on either side */
    @Option(name = "-x", aliases = { "--extend" }, metaVar = "20",
            usage = "distance to extend the genome hit on either side")
    private int extend;

    /** maximum permissible e-value */
    @Option(name = "--maxE", aliases = { "--evalue" }, metaVar = "1e-20", usage = "maximum permissible e-value for a match")
    private double eValue;

    /** temporary directory for BLAST database */
    @Option(name = "--tempDir", metaVar = "Tmp", usage = "temporary directory for BLAST databases")
    private File tempDir;

    /** maximum gap between adjacent sequences */
    @Option(name = "--maxGap", metaVar = "1000", usage = "maximum gap between hits to join")
    private int maxGap;

    /** protein filter file */
    @Option(name = "--filter", metaVar = "proteins.faa", usage = "protein database for filtering input sequences")
    private File filterFile;

    /** protein output file */
    @Option(name = "--fasta", metaVar = "output.faa", usage = "output FASTA file for proteins found")
    private File fastaFile;

    /** upstream region size */
    @Option(name = "--upstream", metaVar = "50", usage = "upstream DNA to include in FASTA output")
    private int upstream;

    /** file containing the target genome */
    @Argument(index = 0, metaVar = "genome.gto", usage = "target genome file containing the proteins",
            required = true)
    private File genomeFile;

    @Override
    protected void setDefaults() {
        this.inFile = null;
        this.batchSize = 10;
        this.eValue = 1e-100;
        this.maxGap = 500;
        this.extend = 50;
        this.tempDir = new File(System.getProperty("user.dir"), "Temp");
        this.filterFile = null;
        this.fastaFile = null;
        this.fastaOut = null;
        this.upstream = 40;
    }

    @Override
    protected boolean validateParms() throws IOException {
        // Read the genome.
        log.info("Loading genome from {}.", this.genomeFile);
        this.genome = new Genome(this.genomeFile);
        // Connect to the DNA input stream.
        if (this.inFile == null) {
            this.inStream = new DnaInputStream(System.in, this.genome.getGeneticCode());
            log.info("DNA sequences will be read from standard input.");
        } else {
            this.inStream = new DnaInputStream(this.inFile, this.genome.getGeneticCode());
            log.info("DNA sequences will be read from {}.", this.inFile);
        }
        if (! this.tempDir.isDirectory()) {
            log.info("Creating temporary file directory {}.", this.tempDir);
            FileUtils.forceMkdir(this.tempDir);
        }
        // Validate the filter file.
        if (this.filterFile == null) {
            log.info("Input sequences will not be filtered.");
            this.filterDb = null;
        } else {
            log.info("Input sequences will be filtered by BLAST database in {}.", this.filterFile);
            try {
                this.filterDb = ProteinBlastDB.createOrLoad(this.filterFile);
            } catch (InterruptedException e) {
                throw new IOException("Interruption in BLAST database creation: " + e.getMessage());
            }
        }
        // Validate the FASTA file.
        if (this.fastaFile != null) {
            log.info("FASTA output sequences will be written to {}.", this.fastaFile);
            this.fastaOut = new FastaOutputStream(this.fastaFile);
        }
        // Validate the parameters.
        if (this.eValue >= 1.0)
            throw new IllegalArgumentException("Invalid eValue specified.  Must be less than 1.");
        if (this.batchSize <= 0)
            throw new IllegalArgumentException("Batch size must be 1 or more.");
        if (this.maxGap < 0)
            throw new IllegalArgumentException("Maximum gap must be 0 or more.");
        if (this.upstream < 0)
            throw new IllegalArgumentException("Upstream length must be 0 or more.");
        if (this.extend < 0)
            throw new IllegalArgumentException("Extension length must be 0 or more.");
        if (! this.genomeFile.canRead())
            throw new FileNotFoundException("Input file" + this.genomeFile + " is not found or unreadable.");
        return true;
    }

    @Override
    public void run() {
        try {
            // Create the BLAST database.
            BlastDB blastDB = this.createBlastDb();
            // Create the BLAST parameters.
            this.mainParms = new BlastParms().maxE(this.eValue);
            this.filterParms = this.mainParms.maxPerQuery(1);
            // Now we loop through the input FASTA stream, building batches.
            int batchCount = 0;
            int seqCount = 0;
            Iterator<SequenceDataStream> batcher = this.inStream.batchIterator(batchSize);
            while (batcher.hasNext()) {
                DnaDataStream batch = (DnaDataStream) batcher.next();
                batchCount++;
                seqCount += batch.size();
                log.info("Processing input batch {} with {} sequences.", batchCount, batch.size());
                this.processBatch(batch, blastDB);
            }
            log.info("All done. {} sequences in {} batches processed.", seqCount, batchCount);
        } catch (Exception e) {
            e.printStackTrace(System.err);
        } finally {
            this.inStream.close();
            if (this.fastaOut != null)
                this.fastaOut.close();
        }
    }

    /**
     * Process a batch of input sequences and write their output.
     *
     * @param batch		batch of DNA sequences to process
     * @param blastDB	blast database to blast against
     *
     * @throws InterruptedException
     * @throws IOException
     */
    private void processBatch(DnaDataStream batch, BlastDB blastDB) throws IOException, InterruptedException {
        DnaDataStream queryStream;
        if (this.filterDb == null) {
            // Here there is no filtering.
            queryStream = batch;
        } else {
            // Here we have to filter using the filter database.
            Map<String, List<BlastHit>> hitMap = BlastHit.sort(batch.blast(this.filterDb, this.filterParms));
            queryStream = new DnaDataStream(this.batchSize, batch.getGeneticCode());
            for (Sequence seq : batch) {
                if (hitMap.containsKey(seq.getLabel())) {
                    queryStream.add(seq);
                }
            }
            log.info("{} sequences left in batch after filtering.", queryStream.size());
        }
        // Perform the BLAST.
        Map<String, List<BlastHit>> hitMap = BlastHit.sort(queryStream.blast(blastDB, this.mainParms));
        // Get a sorted list of the query sequence IDs that produced results.
        ArrayList<String> queryList = new ArrayList<String>(hitMap.keySet());
        queryList.sort(new NaturalSort());
        // Get a map of query sequence IDs to sequences.
        Map<String, String> qMap = queryStream.stream().collect(Collectors.toMap(Sequence::getLabel, Sequence::getSequence));
        // Get a comparator to sort blast hits by subject location.
        Comparator<BlastHit> sortByLoc = new BlastHit.ByLoc(BlastHit.SUBJECT);
        // Get the DNA translator for the query sequences.
        DnaTranslator xlator = new DnaTranslator(this.genome.getGeneticCode());
        // Loop through the query sequences.
        for (String seqId : queryList) {
            // Get the hits for this query.
            List<BlastHit> hits = hitMap.get(seqId);
            hits.sort(sortByLoc);
            // Separate into operon groups.  We want the longest and we will create a single merged location.
            // This will track the current-longest location.
            Location longest = Location.copy(hits.get(0).getSubjectLoc());
            Location current = longest;
            for (int i = 1; i < hits.size(); i++) {
                Location hitLoc = hits.get(i).getSubjectLoc();
                if (hitLoc.distance(current) <= this.maxGap)
                    current.merge(hitLoc);
                else {
                    if (current.getLength() > longest.getLength())
                        longest = current;
                    current = hitLoc;
                }
            }
            if (current.getLength() > longest.getLength())
                longest = current;
            // Expand the location by the extension length on each side.
            int contigLen = this.genome.getContig(longest.getContigId()).length();
            longest.expand(this.extend, this.extend, contigLen);
            // Extract the DNA.
//            String dna = this.genome.getDna(longest);
            // Now we need to get the proteins.
            String targetDna = Contig.reverse(qMap.getOrDefault(seqId, "").toLowerCase());
            System.out.println(seqId);
            if (this.fastaOut != null) {
                List<Sequence> prots = xlator.opTranslate(targetDna, 1, targetDna.length(), seqId, this.upstream);
                for (Sequence prot : prots)
                    fastaOut.write(prot);
            }
            for (int frm = 1; frm <= 3; frm++) {
                String frameProt = xlator.frameTranslate(targetDna, frm);
                System.out.format("    F%d: %s%n", frm, frameProt);
            }
//            List<String> prots = xlator.operonFrom(targetDna);
//            if (prots.size() > 0)
//                System.out.println(dna + "\t" + StringUtils.join(prots, '\t'));
        }
    }

    /**
     * Read in the genome and create the blast database.
     *
     * @throws IOException
     * @throws InterruptedException
     */
    private BlastDB createBlastDb() throws IOException, InterruptedException {
        File tempFile = File.createTempFile("blast", ".fasta", this.tempDir);
        tempFile.deleteOnExit();
        BlastDB retVal = DnaBlastDB.create(tempFile, genome);
        return retVal;
    }

}
