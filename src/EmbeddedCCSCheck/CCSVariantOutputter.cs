using System;
using System.IO;
using System.Text;

using Bio.Variant;
using Bio;

namespace EmbeddedCCSCheck
{
    public class CCSVariantOutputter
    {
        static readonly string header = "Movie,ZMW,Ref,Pos,Length,Type,AtEnd,RefBP,AltBP,IndelType,InHP,Bases,HPLen,HPChar,Read,Direc,ZScore,llseq,llmut";

        public static string GetHeader() {
            return header;
        }

        public static void OutputVariants(Variant v, string movieName, string holeNumber, string[] readName, string[] readDirections, double[] zScores, double[] LL, double[] mutLL, StreamWriter sw) {
            var common = String.Join (",", movieName, holeNumber.ToString (), v.RefName,
                v.StartPosition.ToString (), v.Length.ToString (), v.Type.ToString (), v.AtEndOfAlignment.ToString ());
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
            common += "," + unique + ",";
            lock (sw) {
                for (int i = 0; i < LL.Length; i++) {
                    sw.Write (common);
                    sw.Write (readName[i]+ "," + readDirections[i] + "," + zScores[i] + "," + LL [i] + ", " + mutLL [i] + "\n");

                }
            }
        }
    }
}

