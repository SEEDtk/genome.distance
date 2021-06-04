/**
 *
 */
package org.theseed.genome.signatures;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import org.theseed.genome.Feature;
import org.theseed.proteins.RoleMap;
import org.theseed.utils.ParseFailureException;

/**
 * This signature classifier determines the classes of a feature using the roles.
 *
 * @author Bruce Parrello
 *
 */
public class RoleSignatureClass extends SignatureClass {

    // FIELDS
    /** definition table for useful roles */
    private RoleMap roleMap;

    public RoleSignatureClass(IParms processor) throws IOException, ParseFailureException {
        super(processor);
        File roleFile = processor.getRoleFile();
        if (roleFile == null)
            throw new ParseFailureException("Role definition file is required for signature class ROLE.");
        if (! roleFile.canRead())
            throw new FileNotFoundException("Role definition file " + roleFile + " is not found or unreadable.");
        this.roleMap = RoleMap.load(roleFile);
        log.info("{} roles read from {}.", this.roleMap.size(), roleFile);
    }

    @Override
    protected void addClasses(Set<String> classes, Feature feat) {
        String function = feat.getPegFunction();
        Feature.usefulRoles(this.roleMap, function).stream().forEach(x -> classes.add(x.getId()));
    }

    @Override
    public Map<String, String> getNames(Collection<String> signatures) {
        Map<String, String> retVal = new HashMap<String, String>(EXPECTED_CLASSES);
        for (String signature : signatures)
            retVal.put(signature, this.roleMap.getName(signature));
        return retVal;
    }

}
