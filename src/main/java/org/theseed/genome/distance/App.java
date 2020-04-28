package org.theseed.genome.distance;

import java.util.Arrays;

import org.theseed.sequence.blast.BlastProcessor;
import org.theseed.sequence.blast.MatchProcessor;
import org.theseed.utils.ICommand;

/**
 * These are various commands related to genome and protein distance.
 *
 * build		create or add to a genome minHash database
 * find			query a genome minHash database
 * compare		compare genomes and compute distance based on annotated protein kmers
 * scan			create a role database of the roles in a directory of genomes
 * genomes		compare genomes and compute distance based on DNA kmers
 * mash			compare genomes and compute distance based on minHash
 * width		determine the error rates for different protein sketch widths
 * sketches		convert proteins into signatures (sketches)
 * tune			test different stage and bucket sizes for a protein minHash
 * match		match DNA sequences to proteins in a genome
 * blast		BLAST a sequence source against a BLAST database
 */
public class App
{
    public static void main( String[] args )
    {
        // Get the control parameter.
        String command = args[0];
        String[] newArgs = Arrays.copyOfRange(args, 1, args.length);
        ICommand processor;
        // Parse the parameters.
        switch (command) {
        case "build" :
            processor = new BuildProcessor();
            break;
        case "find" :
            processor = new FindProcessor();
            break;
        case "compare" :
            processor = new DistanceProcessor();
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
        case "width" :
            processor = new WidthProcessor();
            break;
        case "sketches" :
            processor = new SketchProcessor();
            break;
        case "tune" :
            processor = new TuningProcessor();
            break;
        case "match" :
            processor = new MatchProcessor();
            break;
        case "blast" :
            processor = new BlastProcessor();
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
