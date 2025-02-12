/**
 *
 */
package org.theseed.genome.signatures;

import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import org.theseed.genome.Feature;
import org.theseed.p3api.KeyBuffer;
import org.theseed.p3api.P3Connection;
import org.theseed.p3api.P3Connection.Table;

import com.github.cliftonlabs.json_simple.JsonObject;

/**
 * This signature classifier determines the protein class using the global PATRIC protein family.
 *
 * @author Bruce Parrello
 *
 */
public class PgfamSignatureClass extends SignatureClass {

    public PgfamSignatureClass(IParms processor) {
        super(processor);
    }

    @Override
    protected void addClasses(Set<String> classes, Feature feat) {
        String family = feat.getPgfam();
        if (family != null)
            classes.add(family);
    }

    @Override
    public Map<String, String> getNames(Collection<String> signatures) {
        // Connect to PATRIC.
        P3Connection p3 = new P3Connection();
        // Try to find the family names in the PATRIC family table.
        log.info("Retrieving {} family names from PATRIC.", signatures.size());
        Map<String, JsonObject> families = p3.getRecords(Table.FAMILY, signatures, "family_id,family_product");
        log.info("{} records found for the signature families.", families.size());
        Map<String, String> retVal = new HashMap<String, String>(signatures.size());
        for (String signature : signatures) {
            JsonObject record = families.get(signature);
            // There is always the possibility the family is not in PATRIC, because PATRIC is full of errors and omissions.
            // If all we have at the end is an empty string, we create a default name.
            String name = "";
            if (record != null)
                name = KeyBuffer.getString(record, "family_product");
            if (name.isEmpty())
                name = "Missing function " + signature;
            retVal.put(signature, name);
        }
        return retVal;
    }


}
