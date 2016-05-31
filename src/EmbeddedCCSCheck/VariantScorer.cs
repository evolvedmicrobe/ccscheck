using System;
using System.IO;
using System.Runtime.InteropServices;
using System.Runtime.CompilerServices;


using Bio.Variant;

namespace EmbeddedCCSCheck
{
    public static class VariantScorer
    {
        // In C++
        enum MutationType : int
        {
            DELETION = 0,
            INSERTION = 1,
            SUBSTITUTION = 2,
            ANY_INSERTION = 3,
            ANY_SUBSTITUTION = 4
        };
        
        [DllImport ("__Internal", EntryPoint="ScoreVariant")]
        extern static void ScoreVariant (IntPtr ai, int pos, MutationType type, string bases, double[] outputArray); 
        [DllImport ("__Internal", EntryPoint="GetBaseLineLikelihoods")]
        extern static void GetBaseLineLikelihoods (IntPtr ai, double[] outputArray); 
        [MethodImplAttribute(MethodImplOptions.InternalCall)]
        extern static string[] InternalGetReadNames(IntPtr ai);

        /// <summary>
        /// Get the baseline log likelihoods for the current template.
        /// </summary>
        /// <returns>An array of the baseline LL</returns>
        /// <param name="ai">Pointer to the integrator we are using.</param>
        /// <param name="nReads">Number of scores expected, must be set correctly to avoid overflow.</param>
        public static double[] GetBaselineLL(IntPtr ai, int nReads) {
            double[] LL = new double[nReads];
            GetBaseLineLikelihoods (ai, LL);
            return LL;
        }

        /// <summary>
        /// Using the abstract integrator pointed to by ai, score the variant
        /// </summary>
        /// <param name="ai">Pointer to the Integrator used for scoring.</param>
        /// <param name="vi">The variant we are scoring</param>
        /// <param name="nToScore">Number of scores that will be returned. Used to preallocate array.</param>
        public static double[] Score(IntPtr ai, Variant vi, int nToScore) {
            double[] scores = new double[nToScore];
            int pos = vi.StartPosition;
            MutationType type;
            string bases;
            if (vi.Type == VariantType.INDEL) {
                var indel = vi as IndelVariant;
                pos = pos + 1;
                type = indel.InsertionOrDeletion == IndelType.Insertion ? MutationType.INSERTION : MutationType.DELETION;
                bases = indel.InsertedOrDeletedBases;
            } else {
                type = MutationType.SUBSTITUTION;
                var snp = vi as SNPVariant;
                bases = snp.AltBP.ToString ();
            }

            ScoreVariant (ai, pos, type, bases, scores);
            return scores;
        }

        public static string[] GetReadNames(IntPtr ai) {
            return InternalGetReadNames (ai);
        }
    }
}

