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
    public class QVCalibration : CCSReadMetricsOutputter
    {
        const int MAXQV = 94;
        long[] correctCounts = new long[MAXQV];
        long[] incorrectCounts = new long[MAXQV];
        long[] snpErrorCounts = new long[MAXQV];
        long[] deletionErrorCounts = new long[MAXQV];
        long[] insertionErrorCounts = new long[MAXQV];
        string TreatmentName;
        string CoverageLevel;
        public QVCalibration (string dirname):base(dirname, "qv_calibration.csv")
        {
            
        }

        #region implemented abstract members of CCSReadMetricsOutputter

        public override void ConsumeCCSRead (IQualitativeSequence read, PairwiseAlignedSequence aln, List<Variant> variants,
            string treatmentName, string coverageLevel)
        {
            TreatmentName = treatmentName;
            CoverageLevel = coverageLevel;
            if (variants == null) {
                return;
            }
            var refseq = aln.FirstSequence;
            var qseq = aln.SecondSequence as QualitativeSequence;

            foreach (var v in variants) {
                // Ignore lower case bases
                if (!VariantOutputter.VariantValid(v)) {
                    continue;
                }
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
           for (int i = 0; i < refseq.Count; i++) {
                if (refseq [i] == qseq [i] && ((int)qseq[i]) < 95) {
                    correctCounts [qseq.GetQualityScore (i)]++;
                }
            }            
        }
        public override void Finish ()
        {
            SW.WriteLine ("Treatment,Coverage,QV,Correct,Incorrect,SNP,Deletion,Insertion");
            for (int i = 0; i < correctCounts.Length; i++) {
                SW.WriteLine (String.Join (",", TreatmentName,
                    CoverageLevel,
                    i.ToString (), 
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

