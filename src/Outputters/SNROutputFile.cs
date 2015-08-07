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
    public class SNROutputFile : CCSReadMetricsOutputter
    {
        
        public SNROutputFile (string dirName) : base(dirName, "snrs.csv")
        {
            
            SW.WriteLine ("Movie,ZMW,Channel,SNR");
        }

        #region implemented abstract members of CCSReadMetricsOutputter

        public override void ConsumeCCSRead (PacBioCCSRead read, BWAPairwiseAlignment aln, List<Variant> variants)
        {
            
            var prefix = String.Join (",",
                             read.Movie,
                             read.HoleNumber.ToString ());
            SW.WriteLine (prefix + ",A," + read.SnrA);
            SW.WriteLine (prefix + ",C," + read.SnrC);
            SW.WriteLine (prefix + ",G," + read.SnrG);
            SW.WriteLine (prefix + ",T," + read.SnrT);

        }


        #endregion
    }
}

