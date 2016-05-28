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
    public class ZScoreOutputter : CCSReadMetricsOutputter
    {
        public ZScoreOutputter (string dirname) : base(dirname, "zscores.csv")
        {
            SW.WriteLine ("Movie,ZMW,SubRead,ZScore");
        }

        #region implemented abstract members of CCSReadMetricsOutputter

        public override void ConsumeCCSRead (PacBioCCSRead read, BWAPairwiseAlignment aln, List<Variant> variants)
        {
            int i = 1;
            foreach (var v in read.ZScores) {
                var s = String.Join (",", read.Movie, read.HoleNumber, i.ToString (), v.ToString ());
                SW.WriteLine (s);
                i++;
            }            
        }

        #endregion
    }
}

