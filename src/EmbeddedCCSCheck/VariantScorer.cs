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
        extern static void ScoreVariant (IntPtr ai, int Pos, MutationType type, string bases, double[] outputArray); 

        public static unsafe double[] Score(IntPtr ai, Variant vi, int nToScore) {
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

            ScoreVariant (IntPtr.Zero, pos, type, bases, scores);
            Console.WriteLine ("C# Score "  + scores [0]);
            Console.WriteLine ("C# Score " + scores [1]);
            return scores;
        }

    }
}

