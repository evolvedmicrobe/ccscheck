using System;
using System.IO;
using System.Text;
using System.Collections.Generic;
using System.Linq;

using Bio;
using Bio.Variant;
using Bio.BWA;
using Bio.Extensions;
using Bio.IO.SAM;
using Bio.Algorithms.Alignment;

namespace EmbeddedCCSCheck
{

    public static class VariantOutputter
    {
        static StreamWriter SW;
        static BWAPairwiseAligner aligner;
        public static void SetReferenceFastaAndOutputFile(string fastaFile, string outFile) {
            aligner = new BWAPairwiseAligner (fastaFile, false, false);
            SW = new StreamWriter (outFile);
            SW.WriteLine (CCSVariantOutputter.GetHeader ());
        }


        public static void CloseFileStream() {
            SW.Close ();
        }


        public static void AlignAndCallVariants (string seq, IntPtr ai, string movie, long zmw, int nreads)
        {
            if (SW == null || aligner == null) {
                throw new InvalidProgramException ("Tried to align read before calling variants");
            }
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
                    if (variantsInQueryCoordinates.Count != variants.Count) {
                        throw new Exception ("Variants in one direction did not match variants in the other");
                    }

                    double[] baseLL = VariantScorer.GetBaselineLL (ai, nreads);
                    string[] readNames = VariantScorer.GetReadNames (ai).Select(s => s.Split('/').Last()).ToArray();
                    // Now to pair them up

                    for(int i = 0; i < variants.Count; i++) {
                        var refV = variants [i];
                        var queryV = variantsInQueryCoordinates [i];
                        var scores = VariantScorer.Score (ai, queryV, nreads);
                        CCSVariantOutputter.OutputVariants(refV, String.Empty, zmw.ToString(), readNames, baseLL, scores, SW);

                        // Test code not currently used.
#if FALSE
                        #region TEMPORARY_TEST_CODE
                        string newTpl = null;
                        lock(lockobj) {
                            TestResults.Write("Testing " + queryV.StartPosition );
                            TestResults.Write("Mutation at End? " + queryV.AtEndOfAlignment);
                            TestResults.Write("Template Length " + seq.Length);
                            TestResults.WriteLine();
                            TestResults.Flush();



                        if (queryV.Type == VariantType.SNP) {
                            var snp = queryV as SNPVariant;
                            newTpl = VariantScorer.SeqAfterMutation (ai, queryV);
                                TestResults.Write("SNP");
                        } else {
                            var indel = queryV as IndelVariant;
                            if (indel.InsertedOrDeletedBases.Length == 1) {
                                newTpl = VariantScorer.SeqAfterMutation (ai, queryV);
                            }                        
                            TestResults.Write("INDEL");
                            TestResults.Write(indel.InsertionOrDeletion.ToString());
                        }
                        if (newTpl != null) {
                                
                            var orgCount = variants.Count;
                            var mutSeq = new Sequence (DnaAlphabet.Instance, newTpl);
                            result = aligner.AlignRead (mutSeq) as BWAPairwiseAlignment;
                            var newVars = VariantCaller.LeftAlignIndelsAndCallVariants(result);
                            var newCount = newVars.Item2.Count;

                            if (newCount == (orgCount - 1)) {
                                TestResults.WriteLine ("Successfully tested " + zmw + " " + queryV.ToString());
                            } else {
                                TestResults.WriteLine ("Failed on " + zmw + " " + queryV.ToString ());
                                TestResults.WriteLine("Original");
                                Print(aln);
                                TestResults.WriteLine("Mutated");
                                Print(newVars.Item1);
                            }
                        }
                        }
                        #endregion
#endif
                    }
                } 
            }
         
        }

#if FALSE
        static void Print(IPairwiseSequenceAlignment psa) {
            var aln = psa as PairwiseSequenceAlignment;
            TestResults.WriteLine (aln.PairwiseAlignedSequences [0].FirstSequence.ConvertToString ());
            TestResults.WriteLine (aln.PairwiseAlignedSequences [0].SecondSequence.ConvertToString ());
        }
#endif
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
    } 
}

