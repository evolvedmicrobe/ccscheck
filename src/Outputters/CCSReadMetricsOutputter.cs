using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;

using Bio;
using Bio.IO.PacBio;
using Bio.Variant;
using Bio.BWA;
using Bio.BWA.MEM;
using Bio.Algorithms.Alignment;
namespace ccscheck
{
    public abstract class CCSReadMetricsOutputter
    {
        protected StreamWriter SW;
        public CCSReadMetricsOutputter(string dirname, string fname)
        {
            fname = Path.Combine (dirname, fname);
            SW = new StreamWriter (fname);
        }
        public abstract void ConsumeCCSRead(IQualitativeSequence read, PairwiseAlignedSequence aln, List<Variant> variants,
            string treatmentName, string coverageLevel);
        public virtual void Finish() {
            SW.Close ();
        }
    }
}

