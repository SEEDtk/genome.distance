/**
 *
 */
package org.theseed.genome.distance;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.HashSet;
import java.util.Set;

import org.kohsuke.args4j.Argument;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.BaseProcessor;
import org.theseed.counters.CountMap;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.genome.GenomeDirectory;
import org.theseed.proteins.Role;
import org.theseed.proteins.RoleScanner;

/**
 * Scan the genomes in a directory to create a role map.  The role map will be output to a user-specified
 * file and can then be used via the -R option of the DistanceProcessor.
 *
 * The positional parameters are the name of the input genome directory and the name of the output role
 * file.
 *
 * The standard output will contain a report of how often each role occurs.
 *
 * The command-line options are as follows.
 *
 * -v	show more detail on the log
 * -h	display usage information
 *
 *
 * @author Bruce Parrello
 *
 */
public class RoleScanningProcessor extends BaseProcessor {

    // FIELDS

    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(RoleScanningProcessor.class);


    // COMMAND-LINE OPTIONS

    /** input genome directory */
    @Argument(index = 0, metaVar = "genomeDir", usage = "directory of genomes to scan", required = true)
    private File inputDir;

    /** output role file */
    @Argument(index = 1, metaVar = "roleFile", usage = "output role file", required = true)
    private File outFile;

    @Override
    protected void setDefaults() { }

    @Override
    protected boolean validateParms() throws IOException {
        if (! this.inputDir.isDirectory())
            throw new FileNotFoundException("Input directory " + this.inputDir + " is not found or invalid.");
        return true;
    }

    @Override
    public void runCommand() throws Exception {
        RoleScanner roleMap = new RoleScanner();
        GenomeDirectory gDir = new GenomeDirectory(this.inputDir);
        roleMap.addGenomes(gDir);
        log.info("Saving role map to {}.", this.outFile);
        roleMap.save(this.outFile);
        // Now we run through the genomes again, counting each role.
        CountMap<Role> roleCounts = new CountMap<Role>();
        double gCount = 0;
        for (Genome genome : gDir) {
            log.info("Counting roles in {}.", genome);
            // This set will track the unique roles found.  We only want to count a role
            // once per genome.
            Set<Role> roleSet = new HashSet<Role>();
            for (Feature fid : genome.getPegs())
                for (Role role : fid.getUsefulRoles(roleMap))
                    roleSet.add(role);
            // Now count each role.
            for (Role role : roleSet)
                roleCounts.count(role);
            // Count the genome.
            gCount++;
        }
        // Finally, we output the counts.
        log.info("{} roles counted in {} genomes.", roleCounts.size(), gCount);
        System.out.println("Role ID\tRole Name\tCount\tPercent");
        for (CountMap<Role>.Count count : roleCounts.sortedCounts()) {
            Role role = count.getKey();
            double percent = count.getCount() * 100 / gCount;
            System.out.format("%s\t%s\t%6d\t%8.2f%n", role.getId(), role.getName(), count.getCount(),
                    percent);
        }

    }

}
