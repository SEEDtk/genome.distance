package org.theseed.genome.distance;

import java.util.Arrays;

import org.theseed.basic.BaseProcessor;
import org.theseed.genome.signatures.SignatureProcessor;

/**
 * These are various commands related to genome and protein distance.
 *
 * build		create or add to a genome minHash database
 * find			query a genome minHash database
 * scan			create a role database of the roles in a directory of genomes
 * genomes		compare genomes and compute distance based on DNA kmers
 * mash			compare genomes and compute distance based on minHash
 * width		determine the error rates for different protein sketch widths
 * sketches		convert proteins into signatures (sketches)
 * tune			test different stage and bucket sizes for a protein minHash
 * sig			find protein signatures between genome groups
 * methods		compute distances using multiple methods
 * taxCheck		analyze a methods report to determine how each method varies within taxonomic groupings
 * augment		augment a list of pairings to add more taxonomic grouping examples
 * pairs		create a list of pairings from a repgen list file
 * pairMerge	merge pairings from two input files with identical formats
 * methodCorr	compute the correlations between distances in a "methods" report and scores
 * kmerCount	count the number of times each kmer occurs in a set of proteins
 * distCheck	determine how different distance measures relate to groupings
 * distReps		choose representative genomes based on genome distance
 * fastaDist	compute the kmer distance between sequences in a FASTA file
 * fastaReps	compute the number of distance-based representative sequences in a FASTA file
 *
 */
public class App
{
    public static void main( String[] args )
    {
        // Get the control parameter.
        String command = args[0];
        String[] newArgs = Arrays.copyOfRange(args, 1, args.length);
        BaseProcessor processor;
        // Parse the parameters.
        switch (command) {
        case "build" :
            processor = new BuildProcessor();
            break;
        case "find" :
            processor = new FindProcessor();
            break;
        case "scan" :
            processor = new RoleScanningProcessor();
            break;
        case "genomes" :
            processor = new GenomeProcessor();
            break;
        case "mash" :
            processor = new MashProcessor();
            break;
        case "methods" :
            processor = new MethodTableProcessor();
            break;
        case "width" :
            processor = new WidthProcessor();
            break;
        case "sketches" :
            processor = new SketchProcessor();
            break;
        case "tune" :
            processor = new TuningProcessor();
            break;
        case "sig" :
            processor = new SignatureProcessor();
            break;
        case "taxCheck" :
            processor = new TaxCheckProcessor();
            break;
        case "augment" :
            processor = new AugmentProcessor();
            break;
        case "basicPairs" :
            processor = new BasicPairsProcessor();
            break;
        case "pairs" :
            processor = new PairCreateProcessor();
            break;
        case "pairMerge" :
            processor = new PairMergeProcessor();
            break;
        case "kmerCount" :
            processor = new KmerCountProcessor();
            break;
        case "distCheck" :
        	processor = new DistanceCheckProcessor();
        	break;
        case "distReps" :
        	processor = new DistanceRepsProcessor();
        	break;
        case "fastaDist" :
        	processor = new FastaDistanceProcessor();
        	break;
        case "fastaReps" :
        	processor = new FastaDistanceRepsProcessor();
        	break;
        default :
            throw new IllegalArgumentException("Invalid command " + command);
        }

        boolean ok = processor.parseCommand(newArgs);
        if (ok) {
            processor.run();
        }
    }
}
