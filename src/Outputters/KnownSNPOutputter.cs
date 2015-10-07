using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;

using Bio;
using Bio.IO.PacBio;
using Bio.Variant;
using Bio.BWA;
using Bio.BWA.MEM;


namespace ccscheck
{
    public class KnownSNPOutputter : CCSReadMetricsOutputter
    {
        List<KnownSNP> knownSNPs;
        public KnownSNPOutputter (string dirname, BWAPairwiseAligner bwa ):base(dirname, "knownSNPs.csv")
        {
            knownSNPs = KnownVariants.GenerateVariants(bwa).ToList();

        }

        #region implemented abstract members of CCSReadMetricsOutputter

        public override void ConsumeCCSRead (Bio.IO.PacBio.PacBioCCSRead read, Bio.BWA.BWAPairwiseAlignment aln, System.Collections.Generic.List<Bio.Variant.Variant> variants)
        {
            if (aln != null) {
                knownSNPs.ForEach (s => s.ConsumeRead (aln));
            }
        }
        public override void Finish ()
        {
            SW.WriteLine (KnownSNP.ReturnHeader());
            knownSNPs.ForEach (s => SW.WriteLine (s.ReturnDataLine ()));
            base.Finish ();
        }


        #endregion
    }
}

