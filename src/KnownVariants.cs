using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using Bio;
using Bio.IO.PacBio;
using Bio.Variant;
using Bio.BWA;
using Bio.BWA.MEM;
using System;

namespace ccscheck
{
    public class KnownVariants
    {
        public static IEnumerable<KnownSNP> GenerateVariants(BWAPairwiseAligner bwa) {
            for(int i=0; i < variants.Count; i+=2) {
                yield return new KnownSNP (variants[i], variants[i+1], bwa);
            }
        }

        static List<string> variants = new List<string> () {
            // The rest of these are in the file BCRData.txt, can be put online later
            "Q351H",    "AGAAGAAAATCAGCCA[T/C]GGAGTACCTGGAGAAGAAAACACTCAGATCTCG",

        };
    }
}

