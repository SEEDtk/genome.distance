package org.theseed.genome.distance;


/**
 * Hello world!
 *
 */
public class App
{
    public static void main( String[] args )
    {
        DistanceProcessor runObject = new DistanceProcessor();
        boolean ok = runObject.parseCommand(args);
        if (ok) {
            runObject.run();
        }
    }
}
