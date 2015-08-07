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

            SW.WriteLine ("Movie,ZMW,Ref,Pos,Length,Type,QV,AtEnd,RefBP,AltBP,IndelType,InHP,Bases,HPLen,HPChar");

        }

        #region implemented abstract members of CCSReadMetricsOutputter

        public override void ConsumeCCSRead (PacBioCCSRead read, BWAPairwiseAlignment aln, List<Variant> variants)
        {
            if (variants == null) {
                return;
            }
            foreach (var v in variants) {
                var common = String.Join (",", read.Movie, read.HoleNumber.ToString (), v.RefName,
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

