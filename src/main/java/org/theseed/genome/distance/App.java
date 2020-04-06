package org.theseed.genome.distance;

import java.util.Arrays;

import org.theseed.utils.ICommand;

/**
 * Hello world!
 *
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
        case "compare" :
            processor = new DistanceProcessor();
            break;
        case "scan" :
            processor = new RoleScanningProcessor();
            break;
        case "genomes" :
            processor = new GenomeProcessor();
            break;
        case "sketches" :
            processor = new MashProcessor();
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
