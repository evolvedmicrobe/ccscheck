using System;
using System.IO;
using System.Text;
using System.Collections.Generic;
using System.Collections.Concurrent;
using System.Threading.Tasks;
using System.Linq;

using Bio;
using Bio.IO;
using Bio.IO.PacBio;
using Bio.Variant;
using Bio.BWA.MEM;
using Bio.BWA;
using Bio.Extensions;
using Bio.IO.SAM;
using Bio.Algorithms.Alignment;

namespace EmbeddedCCSCheck
{
    public enum EventType : int { SNP=0, Insertion=1, Deletion=2};

    public class VariantInfo {
        public string RefVariantDescription;
        public EventType Type;
        public Variant TemplateMutation;
        public VariantInfo(Variant refCall, Variant queryCall) {            
            // Load up all the important information about the reference variant.
            var refVariantDescription = new StringBuilder (20);
            refVariantDescription.Append (refCall.RefName);
            refVariantDescription.Append (':');
            refVariantDescription.Append (refCall.StartPosition);
            refVariantDescription.Append (':');
            if (refCall.Type == VariantType.SNP) {
                var snp = refCall as SNPVariant;
                refVariantDescription.Append ("S:" + snp.AltBP.ToString());
            } else if (refCall.Type == VariantType.INDEL) {
                var indel = refCall as IndelVariant;
                if (indel.InsertionOrDeletion == IndelType.Insertion) {
                    refVariantDescription.Append ("I" + indel.InsertedOrDeletedBases.Length + ":" + indel.InsertedOrDeletedBases);
                } else if (indel.InsertionOrDeletion == IndelType.Deletion) {
                    refVariantDescription.Append ("D" + indel.InsertedOrDeletedBases.Length + ":" + indel.InsertedOrDeletedBases);
                } else {
                    throw new Exception ("Complex types not allowed");
                }
            } else {
                throw new Exception ("Complex variants not allowed");
            }
            this.RefVariantDescription = refVariantDescription.ToString ();

            // Now the ever important variant information
            this.TemplateMutation = queryCall;
        }

      
    }


    public static class Aligner
    {
        public static VariantInfo CreateVariantInfoFromVariant(Variant v, int qStart, bool orgRevComp, SAMAlignedSequence seq) {

            return null;
        }
        static BWAPairwiseAligner aligner;

        public static void SetReferenceFasta(string fastaFile) {
            aligner = new BWAPairwiseAligner (fastaFile, false, false);
        }
        public static void Align (string seq)
        {
            var read = new Sequence(DnaAlphabet.Instance, seq);
            var result = aligner.AlignRead (read) as BWAPairwiseAlignment;
            if (result == null) {
                return;
            } else {
                var alnv = VariantCaller.LeftAlignIndelsAndCallVariants(result);
                var variants = alnv.Item2;
                var aln = alnv.Item1;
                // Super annoying, I need to get back the query positions of this stuff 
                // to test everything, which is information I specifically discarded before
                if (variants.Count > 0) {
                    // TODO: Super lazy here, we'll get the variants in query coordinates by flipping,
                    // and recalling everything, this is a very inefficient operation.
                    var pairwise = new PairwiseSequenceAlignment ();
                    var seqpairs = new PairwiseAlignedSequence ();
                    if (result.OriginallyReverseComplemented) {
                        variants.Reverse ();
                        seqpairs.FirstSequence = result.AlignedQuerySeq.GetReverseComplementedSequence ();
                        seqpairs.SecondSequence = result.AlignedRefSeq.GetReverseComplementedSequence ();
                    } else {
                        seqpairs.FirstSequence = result.AlignedQuerySeq;
                        seqpairs.SecondSequence = result.AlignedRefSeq;
                    }
                    pairwise.Add (seqpairs);
                    var variantsInQueryCoordinates = VariantCaller.CallVariants(pairwise);
                    var se = GetQueryStartAndEndPadding(result.AlignedSAMSequence);
                    var offset = result.OriginallyReverseComplemented ? se.Item2 : se.Item1;
                    variantsInQueryCoordinates.ForEach (v => v.StartPosition += offset);
                    variants.ForEach (p => {
                        p.StartPosition += result.AlignedSAMSequence.Pos;
                        p.RefName = result.Reference;
                    });
                    var results = Enumerable.Zip (variants, variantsInQueryCoordinates, (x, y) => new VariantInfo (x, y)).ToList(); 

                    foreach (var v in variantsInQueryCoordinates) {
                        VariantScorer.Score (IntPtr.Zero, v, 10);
                    }
                } 
            }
        }
        static Tuple<int, int> GetQueryStartAndEndPadding(SAMAlignedSequence seq) {
           
            var cigar = CigarUtils.GetCigarElements (seq.CIGAR);
            int start = GetFirstNonClippedPosition (cigar);
            cigar.Reverse ();
            int end = GetFirstNonClippedPosition (cigar);
            return new Tuple<int, int> (start, end);
        }

        static int GetFirstNonClippedPosition(IEnumerable<CigarElement> cigar) {
            int start = 0;
            foreach (var c in cigar) {
                char type = c.Operation;
                if (type == CigarOperations.PADDING ||
                    type == CigarOperations.HARD_CLIP ||
                    type == CigarOperations.SKIPPED) {
                    throw new FormatException ("Padding, hard clipping and skipped bases are not currently supported.");
                } else if (type == CigarOperations.SOFT_CLIP) {
                    start += c.Length;
                } else {
                    return start;
                }
            }
            return start;

        }

        
    } 
}

