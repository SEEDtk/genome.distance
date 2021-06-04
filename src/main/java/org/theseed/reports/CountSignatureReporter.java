/**
 *
 */
package org.theseed.reports;

import java.io.PrintWriter;

/**
 * This is the simplest signature report:  it shows the signature classes along with the number and percent each
 * occurs in the genomes.
 *
 * @author Bruce Parrello
 *
 */
public class CountSignatureReporter extends SignatureReporter {

    public CountSignatureReporter(PrintWriter writer, IParms processor) {
        super(writer, processor);
    }

    @Override
    public void initReport() {
        // Write the header line.
        this.writer.println("class\tclass_name\tin_count\tout_count\tin_percent\tout_percent");
    }

    @Override
    public void showClass(String signature, int count1, int count2) {
        this.writer.format("%s\t%s\t%d\t%d\t%5.1f\t%5.1f%n", signature, this.getName(signature), count1, count2,
                this.getPercent1(count1), this.getPercent2(count2));
    }

    @Override
    public void closeReport() {
    }

}
