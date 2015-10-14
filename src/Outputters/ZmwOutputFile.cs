using System.IO;
using System.Collections.Generic;
using System.Linq;
using System;

using Bio;
using Bio.IO.PacBio;
using Bio.Variant;
using Bio.BWA;
using Bio.BWA.MEM;

namespace ccscheck
{
    public class ZmwOutputFile : CCSReadMetricsOutputter
    {
        public ZmwOutputFile (string dirname) : base(dirname, "zmws.csv")
        {
            SW.WriteLine ("Movie,ZMW,SnrA,SnrC,SnrG,SnrT,Aligned,RQ,AvgZScore,ReadLength,NP,NumSuccess,NumZlow,NumABMis,NumOther,AlnLength,NumErrors,NumIndel,NumSNP,Ref");
        }

        #region implemented abstract members of CCSReadMetricsOutputter

        public override void ConsumeCCSRead (PacBioCCSRead read, BWAPairwiseAlignment aln, List<Variant> variants)
        {
            var aligned = aln==null ? "FALSE" : "TRUE";
            string start = String.Join (",", read.Movie, read.HoleNumber.ToString (),
                               read.SnrA.ToString (), read.SnrC.ToString (), read.SnrG.ToString (),
                               read.SnrT.ToString (), aligned, read.ReadQuality.ToString (),
                               read.AvgZscore.ToString (), read.Sequence.Count.ToString (), read.NumPasses.ToString (),
                               read.ReadCountSuccessfullyAdded.ToString(), read.ReadCountBadZscore.ToString(), 
                               read.ReadCountAlphaBetaMismatch.ToString(),
                read.ReadCountOther.ToString());
            if (aln == null) {
                start = start + ",NA,NA,NA,NA,NA";
            } else {
                var indels = variants.Count (z => z.Type == VariantType.INDEL);
                var snps = variants.Count (z => z.Type == VariantType.SNP);
                var total = variants.Count;
                var length = aln.AlignedSAMSequence.RefEndPos - aln.AlignedSAMSequence.Pos;
                start = start + "," + String.Join (",", length.ToString (), total.ToString (), indels.ToString (), snps.ToString (),
                    aln.AlignedSAMSequence.RName);
            }
            SW.WriteLine (start);
        }

       
        #endregion
    }
}

