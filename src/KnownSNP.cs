using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Bio;
using Bio.IO.PacBio;
using Bio.Variant;
using Bio.BWA;
using Bio.BWA.MEM;
using System.Threading;

namespace ccscheck
{
    public class KnownSNP
    {
        public const int MIN_QUAL_VALUE = 45;
        public readonly string Name;
        public readonly string VariantString;
        public readonly char WT;
        public readonly char ALT;
        readonly char LeftBase;
        readonly char RightBase;
        int start, end;
        string ChrName;
        int snpPosition;
        public int ReadsSeen, ReadsOverlap, ReadsWT, ReadsAlt, ReadsOther;
        /// <summary>
        /// Initializes a new instance of the <see cref="ccscheck.CMLVariant"/> class.
        ///  "M351T",    "CACTCAGATCTCGTCAGCCA[T/C]GGAGTACCTGGAGAAGAAAA",
        /// </summary>
        /// <param name="name">Name.</param>
        /// <param name="sequence">Sequence.</param>
        public KnownSNP (string name, string sequence, BWAPairwiseAligner bwa)
        {
            Name = name;
            VariantString = sequence;
            var leftStart = sequence.IndexOf ('[');
            var rightStart = sequence.IndexOf (']');
            WT = sequence [leftStart + 1];
            ALT = sequence [rightStart - 1];
            LeftBase = sequence [leftStart - 1];
            RightBase = sequence [rightStart + 1];
            var fullSeq = sequence.Substring (0, leftStart) + WT + sequence.Substring (rightStart + 1, sequence.Length - rightStart -1); 
            var seq = new Sequence (DnaAlphabet.Instance, fullSeq);
            var aln = bwa.AlignRead (seq) as BWAPairwiseAlignment;
            if (aln == null) {
                throw new Exception ("Variant: " + name + " - " + sequence + " did not align to the reference."); 
            }
            if (aln.AlignedSAMSequence.CIGAR != (seq.Count.ToString () + "M")) {
                throw new Exception ("Variant: " + name + " - " + sequence + " was not a perfect match to the reference.  Cigar was: " + aln.AlignedSAMSequence.CIGAR);
            }
            if (aln.OriginallyReverseComplemented) {
                throw new Exception ("Variant aligns as a reverse complement relative to reference");
            }
            ChrName = aln.AlignedSAMSequence.RName;
            start = aln.AlignedSAMSequence.Pos;
            end = aln.AlignedSAMSequence.RefEndPos;
            snpPosition = start + leftStart; 
        }

        /// <summary>
        /// Consumes an alignment and adds it to counting statistics in a thread safe fashion.
        /// </summary>
        /// <param name="aln">Aln.</param>
        public void ConsumeRead(BWAPairwiseAlignment aln) {
            if(aln != null) {
                Interlocked.Increment(ref ReadsSeen);
                // Verify overlaps with region
                if (aln.Reference ==  this.ChrName &&
                    aln.AlignedSAMSequence.Pos < start &&
                    aln.AlignedSAMSequence.RefEndPos > end) 
                {
                    Interlocked.Increment(ref ReadsOverlap);
                    // Now go find SNP
                    var ref_pos = aln.AlignedSAMSequence.Pos - 1;
                    var refSeq = aln.AlignedRefSeq;
                    var querySeq = aln.AlignedQuerySeq as QualitativeSequence;
                    int i = 0;
                    bool found = false;
                    for (; i < refSeq.Count; i++) {
                        if (refSeq [i] != '-') {
                            ref_pos++;
                            if (ref_pos == snpPosition) {
                                found = true;
                                break;
                            }
                        }
                    }
                    if (!found) {
                        throw new Exception ("Could not find SNP in overlapping alignment");
                    }
//                    Console.WriteLine (i);
//                    Console.WriteLine (refSeq.ConvertToString ());
//                    Console.WriteLine ((new Sequence(querySeq)).ConvertToString ());
//                    Console.WriteLine (refSeq [i - 1]);
//                    Console.WriteLine (refSeq [i]);
//                    Console.WriteLine (refSeq [i + 1]);
//                    Console.WriteLine(querySeq.GetPhredQualityScore(i - 1));
//                    Console.WriteLine(querySeq.GetPhredQualityScore(i));
//                    Console.WriteLine(querySeq.GetPhredQualityScore(i + 1));
                    // Check the 2 neighboring bases match and all 3 have QV values above the minimum
                    if (refSeq [i - 1] == LeftBase &&
                        refSeq [i] == WT &&
                        refSeq [i + 1] == RightBase &&                       
                        querySeq[i - 1] == LeftBase && 
                        querySeq[i + 1] == RightBase && 
                        querySeq.GetPhredQualityScore(i - 1) >= MIN_QUAL_VALUE &&
                        querySeq.GetPhredQualityScore(i + 1) >= MIN_QUAL_VALUE &&
                        querySeq.GetPhredQualityScore(i) >= MIN_QUAL_VALUE)
                    {
                        var bp = querySeq [i];
                        if (bp == WT) {
                            Interlocked.Increment (ref ReadsWT);
                        } else if (bp == ALT) {
                            Console.WriteLine (i);
                            Console.WriteLine (refSeq.ConvertToString ());
                            Console.WriteLine ((new Sequence(querySeq)).ConvertToString ());
                            Console.WriteLine (refSeq [i - 1]);
                            Console.WriteLine (refSeq [i]);
                            Console.WriteLine (refSeq [i + 1]);
                            Console.WriteLine(querySeq.GetPhredQualityScore(i - 1));
                            Console.WriteLine(querySeq.GetPhredQualityScore(i));
                            Console.WriteLine(querySeq.GetPhredQualityScore(i + 1));

                            Interlocked.Increment (ref ReadsAlt);
                        } else if (bp == ReadsOther) {
                            Interlocked.Increment (ref ReadsOther);
                        }
                    }
                }
            }
        }
        public static string ReturnHeader() {
            return "Name,Seq,Start,Total,Overlap,WT,Alt,Other,PercentAlt";
        }
        public string ReturnDataLine() {
            return String.Join (",", Name, VariantString, snpPosition.ToString (), ReadsSeen.ToString (),
                ReadsOverlap.ToString (), ReadsWT.ToString (), ReadsAlt.ToString (),
                ReadsOther.ToString (), (ReadsAlt / (double)(ReadsAlt + ReadsWT + ReadsOther)).ToString ());
        }
    }
}

