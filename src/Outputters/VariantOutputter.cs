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
    public class VariantOutputter : CCSReadMetricsOutputter
    {
        public VariantOutputter (string dirname) : base (dirname, "variants.csv")
        {

            SW.WriteLine ("Treatment,Coverage,Ref,Pos,Length,Type,QV,AtEnd,RefBP,AltBP,IndelType,InHP,Bases,HPLen,HPChar");

        }

        #region implemented abstract members of CCSReadMetricsOutputter
        /// <summary>
        /// We should never call variants where the reference is different, but this checks for insertions and snps 
        /// that we don't, unfortunately we can never check Deletions without more code modifications.
        /// </summary>
        /// <returns><c>true</c>, if valid was varianted, <c>false</c> otherwise.</returns>
        /// <param name="v">V.</param>
        public static bool VariantValid(Variant v) {
            var snp = v as SNPVariant;
            var indel = v as IndelVariant;
            if (snp != null) {                
                if (((int)snp.AltBP) > 95) {
                    throw new Exception("Error should not be called in uncovered areas!");
                }
            } else if (indel != null) {
                if(indel.InsertionOrDeletion == IndelType.Insertion) {
                    if (((int)indel.InsertedOrDeletedBases [0]) > 95) {
                        throw new Exception("Error should not be called in uncovered areas!");
                    }
                }               
            } else {
                throw new Exception ("Variant is neither SNP nor Indel");
            }
            return true;
        }
        public override void ConsumeCCSRead (IQualitativeSequence read, BWAPairwiseAlignment aln, List<Variant> variants,
            string treatmentName, string coverageLevel)
        {
            if (variants == null) {
                return;
            }
            var refseq = aln.AlignedRefSeq;

            foreach (var v in variants) {
                // Ignore lowercase bases
                if (!VariantValid(v)) {
                    continue;
                }
                var common = String.Join (",", treatmentName, coverageLevel, v.RefName,
                                 v.StartPosition.ToString (), v.Length.ToString (), v.Type.ToString (), v.QV.ToString (), v.AtEndOfAlignment.ToString ());
                var snp = v as SNPVariant;
                string unique;
                if (snp != null) {
                    unique = String.Join (",", snp.RefBP, snp.AltBP, "NA,NA,NA,NA,NA");
                } else {
                    var indel = v as IndelVariant;
                    unique = string.Join (",", "NA,NA", indel.InsertionOrDeletion.ToString (),
                        indel.InHomopolymer.ToString (), indel.InsertedOrDeletedBases.ToString (),
                        indel.HomopolymerLengthInReference.ToString (), indel.HomopolymerBase);
                }
                SW.WriteLine (common + "," + unique);
            }
        }

        #endregion
    }
}

