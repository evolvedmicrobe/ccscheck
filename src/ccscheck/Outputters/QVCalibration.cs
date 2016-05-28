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
    public class QVCalibration : CCSReadMetricsOutputter
    {
        const int MAXQV = 94;
        long[] correctCounts = new long[MAXQV];
        long[] incorrectCounts = new long[MAXQV];
        long[] snpErrorCounts = new long[MAXQV];
        long[] deletionErrorCounts = new long[MAXQV];
        long[] insertionErrorCounts = new long[MAXQV];

        public QVCalibration (string dirname):base(dirname, "qv_calibration.csv")
        {
            
        }

        #region implemented abstract members of CCSReadMetricsOutputter

        public override void ConsumeCCSRead (PacBioCCSRead read, BWAPairwiseAlignment aln, List<Variant> variants)
        {
            if (variants == null) {
                return;
            }
            foreach (var v in variants) {                
                incorrectCounts [v.QV]++;
                if (v.Type == VariantType.SNP) {
                    snpErrorCounts [v.QV]++;

                } else {
                    var id = v as IndelVariant;
                    if (id.InsertionOrDeletion == IndelType.Deletion) {
                        deletionErrorCounts [v.QV]++;
                    } else if (id.InsertionOrDeletion == IndelType.Insertion) {
                        insertionErrorCounts [v.QV]++;
                    }
                }
            }
            var refseq = aln.AlignedRefSeq;
            var qseq = aln.AlignedQuerySeq as QualitativeSequence;
            for (int i = 0; i < refseq.Count; i++) {
                if (refseq [i] == qseq [i]) {
                    correctCounts [qseq.GetQualityScore (i)]++;
                }
            }            
        }
        public override void Finish ()
        {
            SW.WriteLine ("QV,Correct,Incorrect,SNP,Deletion,Insertion");
            for (int i = 0; i < correctCounts.Length; i++) {
                SW.WriteLine (String.Join (",", i.ToString (), 
                    correctCounts [i].ToString (), 
                    incorrectCounts [i].ToString (),
                    snpErrorCounts[i].ToString(),
                    deletionErrorCounts[i].ToString(),
                    insertionErrorCounts[i].ToString()));
            }
            base.Finish ();
        }

        #endregion
    }
}

