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
    public class AlnPair {
        public PacBioCCSRead Read;
        public BWAPairwiseAlignment Aln;
        public AlnPair(PacBioCCSRead read, BWAPairwiseAlignment aln)
        {
            this.Read = read;
            this.Aln = aln;
        }
    }
    public class GenotypeOutputter : CCSReadMetricsOutputter
    {
        int brca_start = 32889619;
        int brca_end = 32973810;
        List<AlnPair> alns = new List<AlnPair> ();

        public GenotypeOutputter(string dirname, string fname) : base(dirname, fname)
        {
            
        }

        #region implemented abstract members of CCSReadMetricsOutputter

        public override void ConsumeCCSRead (Bio.IO.PacBio.PacBioCCSRead read, Bio.BWA.BWAPairwiseAlignment aln, System.Collections.Generic.List<Bio.Variant.Variant> variants)
        {
            throw new NotImplementedException ();
        }

        #endregion


    }
}

