using System;
using System.IO;
using System.Runtime.InteropServices;
using System.Runtime.CompilerServices;


using Bio.Variant;

namespace EmbeddedCCSCheck
{
    public class VariantScorer
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

        #region STATICS_AND_CALLERS
        // TODO: Change P/Invoke to internal calls, which are a heck of a lot faster.


        /// Using the abstract integrator pointed to by ai, score the variant
        /// </summary>
        /// <param name="ai">Pointer to the Integrator used for scoring.</param>
        /// <param name="vi">The variant we are scoring</param>
        /// <param name="nToScore">Number of scores that will be returned. Used to preallocate array.</param>
        public double [] Score (Variant vi)
        {
            double [] scores = new double [numEvaluators];
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

            bool noError = ScoreVariant (abstractIntegrator, pos, type, bases, scores);
            if (!noError) {
                throw new Exception ("Failed to test mutation - Exception in C++");
            }
            return scores;
        }
        [DllImport ("__Internal", EntryPoint="ScoreVariant")]
        extern static bool ScoreVariant (IntPtr ai, int pos, MutationType type, string bases, double[] outputArray);

        /// <summary>
        /// Get the baseline log likelihoods for the current template.
        /// </summary>
        /// <returns>An array of the baseline LL</returns>
        public double [] GetBaselineLL ()
        {
            double [] LL = new double [numEvaluators];
            GetBaseLineLikelihoods (abstractIntegrator, LL);
            return LL;
        }
        [DllImport ("__Internal", EntryPoint = "GetBaseLineLikelihoods")]
        extern static void GetBaseLineLikelihoods (IntPtr ai, double [] outputArray);

        /// <summary>
        /// Get the ZScores
        /// </summary>
        /// <returns>An array of the baseline LL</returns>
        public double [] GetZScores ()
        {
            double [] zs = new double [numEvaluators];
            GetBaseLineZScores (abstractIntegrator, zs);
            return zs;
        }
        [DllImport ("__Internal", EntryPoint = "GetBaseLineZScores")]
        extern static void GetBaseLineZScores (IntPtr ai, double [] outputArray);

        public string [] GetReadNames ()
        {
            return InternalGetReadNames (abstractIntegrator);
        }
        [MethodImplAttribute(MethodImplOptions.InternalCall)]
        extern static string[] InternalGetReadNames(IntPtr ai);

        public string [] GetReadDirections ()
        {
            return InternalGetReadDirections (abstractIntegrator);
        }
        [MethodImplAttribute (MethodImplOptions.InternalCall)]
        extern static string [] InternalGetReadDirections (IntPtr ai);

        // Used to verify the application of the mutation removes the variant
        [MethodImplAttribute(MethodImplOptions.InternalCall)]
        extern static string GetTemplateAfterMutation(IntPtr ai, int pos, int type, byte mutBase);

        [MethodImplAttribute (MethodImplOptions.InternalCall)]
        extern static int InternalGetEvaluatorCount (IntPtr ai);

        #endregion

        private IntPtr abstractIntegrator;
        private int numEvaluators;

        public VariantScorer (IntPtr ai)
        {
            abstractIntegrator = ai;
            numEvaluators = InternalGetEvaluatorCount (abstractIntegrator);
        }
        /// <summary>


// Temporary verification code, to be restructured for testing later.
#if FALSE
        public static string SeqAfterMutation(IntPtr ai, Variant vi) {
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
            return GetTemplateAfterMutation (ai, pos, (int)type, Convert.ToByte(bases[0]));
        }
#endif
    }
}

