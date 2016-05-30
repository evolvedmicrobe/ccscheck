using System;
using System.IO;
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

namespace EmbeddedCCSCheck
{
    class MainClass
    {
        public static void Main (string[] args)
        {
            PlatformManager.Services.MaxSequenceSize = int.MaxValue;
            PlatformManager.Services.DefaultBufferSize = 4096;
            PlatformManager.Services.Is64BitProcessType = true;
            //EmbeddedCCSCheck.Aligner.SetReferenceFastaAndOutputFile ("/Users/nigel/BroadBundle/human_g1k_v37.fasta", "tmp.txt");
            //EmbeddedCCSCheck.Aligner.Align("gatgggCCAaccttgtggaagTagaggtgccagCAaatagtctgggaagggga", IntPtr.Zero, "ASDFA" , 5464L, 25);
        }
    }
}
