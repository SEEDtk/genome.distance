package org.theseed.genome.distance;

import java.util.Arrays;

import org.theseed.genome.signatures.SignatureProcessor;
import org.theseed.utils.BaseProcessor;

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
        default :
            throw new IllegalArgumentException("Invalid command " + command);
        }

        boolean ok = processor.parseCommand(newArgs);
        if (ok) {
            processor.run();
        }
    }
}
