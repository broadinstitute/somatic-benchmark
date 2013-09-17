/*
* Copyright (c) 2012 The Broad Institute
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.cga.tools.gatk.walkers.varianteval.stratifications;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.varianteval.stratifications.VariantStratifier;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.exceptions.UserException.BadInput;
import org.broadinstitute.variant.variantcontext.Genotype;
import org.broadinstitute.variant.variantcontext.GenotypesContext;
import org.broadinstitute.variant.variantcontext.VariantContext;

import java.util.Collections;
import java.util.List;

/**
 * Stratifies the eval RODs by the alternate allele read fraction across samples
 *  (calculated as AD of alternate allele / DP)
 *
 *  WARNING: This is intended for a very narrow purpose and combines read fractions across all given genotypes.
 *  WARNING: It explicitly looks at the first alternate allele only;
 *
 * Uses a constant 0.005 frequency grid, and projects the AF INFO field value.  Requires
 * that AF be present in every ROD, otherwise this stratification throws an exception
 */
public class AltReadFraction extends VariantStratifier {
    @Override
    public void initialize() {
        for( double a = 0.000; a <= 1.005; a += 0.005 ) {
            states.add(String.format("%.3f", a));
        }
    }

    public List<Object> getRelevantStates(ReferenceContext ref, RefMetaDataTracker tracker, VariantContext comp, String compName, VariantContext eval, String evalName, String sampleName) {
        if (eval != null) {
            try {

                GenotypesContext genotypes = eval.getGenotypes();
                if (genotypes==null)
                    throw new BadInput("Stratification by AltReadFraction requires genotype information in the vcf");

                int combinedAlternateDepth = 0;
                int combinedDP = 0;

                for(Genotype genotype : genotypes){
                    combinedAlternateDepth += genotype.getAD()[1];
                    combinedDP += genotype.getDP();
                }

                double alternateReadFraction = (double)combinedAlternateDepth / (double)combinedDP;

                return Collections.singletonList((Object)String.format("%.3f", 5.0 * MathUtils.round(alternateReadFraction/5.0, 3)));
            } catch (Exception e) {
                return Collections.emptyList();
            }
        }

        return Collections.emptyList();
    }
}
