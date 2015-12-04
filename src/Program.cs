using System;
using System.IO;
using System.Collections.Generic;
using System.Collections.Concurrent;
using System.Threading.Tasks;
using System.Linq;
using System.Threading;

using Bio;
using Bio.IO.PacBio;
using Bio.Variant;
using Bio.BWA.MEM;
using Bio.IO.FastQ;
using Bio.BWA;
using Bio.Extensions;
using Bio.Algorithms.MUMmer;
using Bio.Algorithms.Alignment;


namespace ccscheck
{
    class MainClass
    {
        public static void Main (string[] args)
        {
            try {
                PlatformManager.Services.MaxSequenceSize = int.MaxValue;
                PlatformManager.Services.DefaultBufferSize = 4096;
                PlatformManager.Services.Is64BitProcessType = true;

                if (args.Length > 3) {
                    Console.WriteLine ("Too many arguments");
                    DisplayHelp ();
                } else if (args.Length < 3) {
                    Console.WriteLine("Not enough arguments");
                    DisplayHelp();
                }else if (args [0] == "h" || args [0] == "help" || args [0] == "?" || args [0] == "-h") {
                    DisplayHelp ();
                } else {
                    string fastq_name = args [0];
                    string out_dir = args [1];
                    string ref_name = args.Length > 2 ? args [2] : null;

                    if (!File.Exists(fastq_name)) {
                        Console.WriteLine ("Can't find file: " + fastq_name);
                        return;
                    }
                    if (ref_name != null && !File.Exists (ref_name)) {
                        Console.WriteLine ("Can't find file: " + ref_name);
                        return;
                    }
                    if (Directory.Exists (out_dir)) {
                        Console.WriteLine ("The output directory already exists, please specify a new directory or delete the old one.");
                        //return;
                    }

                    // Get the Coverage level
                    var cov_level = fastq_name.Split('_').Last().Split('.')[0];
                    var treatment = Path.GetDirectoryName(fastq_name).Split('/').Last();

                    Directory.CreateDirectory (out_dir);

                    List<CCSReadMetricsOutputter> outputters = new List<CCSReadMetricsOutputter> () { 
                        new VariantOutputter(out_dir),
                        new QVCalibration(out_dir)};

                    FastQParser fqp = new FastQParser();
                    fqp.Alphabet = DnaAlphabet.Instance;
                    fqp.FormatType = FastQFormatType.Sanger;
                    MUMmerAligner mum = new MUMmerAligner();

                    var fp = new Bio.IO.FastA.FastAParser();
                    var refseq = fp.Parse(ref_name).First();

                    // Produce aligned reads with variants called in parallel.
                    var reads = new BlockingCollection<Tuple<IQualitativeSequence, PairwiseAlignedSequence, List<Variant>>>();
                    Task producer = Task.Factory.StartNew(() =>
                        {
                            try 
                            {
                                Parallel.ForEach(fqp.Parse(fastq_name), z => {
                                    try {
                                        List<Variant> variants = null;
                                        // Save the lower case version, make a new upper case version for aligning and variant calling
                                        var z2 = new Sequence(z.Alphabet, z.ConvertToString().ToUpper(), true);
                                        z2.ID = z.ID;
                                        var aln = mum.AlignSimple(refseq, new List<Sequence>() { z2}).First();
                                        if (aln != null) {

                                            // Now to add back in the quality scores and lowercase bases
                                            var zasq = z as QualitativeSequence;
                                            var old_seq = aln.PairwiseAlignedSequences[0].SecondSequence as Sequence;
                                            var oss = old_seq.ConvertToString();
                                            var cnt = oss.Replace("-","").Length;
                                            if (cnt != zasq.Count) {
                                               throw new Exception("Super evil!");
                                            }
                                            var new_ref = aln.PairwiseAlignedSequences[0].FirstSequence;                                                
                                            byte[] qual_scores = new byte[old_seq.Count];
                                            byte[] data = new byte[old_seq.Count];
                                            int readPos =0;
                                            for(int i = 0; i < oss.Length; i++)
                                            {
                                                byte newBase, qv;
                                                var cur = oss[i];
                                                if(cur != '-') {
                                                    qv = (byte)zasq.GetQualityScore(readPos);
                                                    newBase = zasq[readPos];
                                                    readPos++;
                                                    if (newBase > 95)
                                                    {
                                                        if (cur != new_ref[i]) {
                                                            throw new Exception("Danger at the boundary!");
                                                        }
                                                    }
                                                    if (cur != newBase && cur != (newBase - 32)) {
                                                        for(int j = i-10; j < i+10; j++) 
                                                        {
                                                            Console.Write((char)oss[j]);
                                                        }
                                                        Console.Write("\n");
                                                        for(int j = readPos-10; j < readPos+10; j++) 
                                                        {
                                                            Console.Write((char)zasq[j]);
                                                        }                                                             
                                                        throw new Exception("Indexing error, new and old are off!");
                                                    }
                                                }
                                                else {
                                                    qv = 0;
                                                    newBase = (byte)cur;
                                                }
                                                qual_scores[i] = (byte)(qv + 33);
                                                data[i] = newBase;
                                            }
                                            var Seq2 = new QualitativeSequence(DnaAlphabet.Instance, zasq.FormatType, old_seq.ToArray(), qual_scores, false);
                                            Seq2.ID = zasq.ID;
                                            aln.PairwiseAlignedSequences[0].SecondSequence = Seq2;
                                            // Finish converting back



                                            variants = VariantCaller.CallVariants(aln);

                                            variants.ForEach( p => {
                                                //p.StartPosition += aln.AlignedSAMSequence.Pos;
                                                p.RefName = aln.FirstSequence.ID;
                                                });
                                        }
                                        var res = new Tuple<IQualitativeSequence, PairwiseAlignedSequence, List<Variant>>(z, aln.PairwiseAlignedSequences.First(), variants);
                                        reads.Add(res);
                                    }
                                    catch(Exception thrown) {
                                        Console.WriteLine("READ FAILED: " + fastq_name +"\nREAD=" + z.ID);
                                        Console.WriteLine(thrown.Message);
                                    } });
                            } catch(Exception thrown) {
                                Console.WriteLine(thrown.StackTrace);
                                Console.WriteLine("Could not parse FASTQ file: " + fastq_name + "\n" + thrown.Message);
                                while(thrown.InnerException!=null) {
                                    Console.WriteLine(thrown.InnerException.StackTrace);
                                    Console.WriteLine(thrown.InnerException.Message);
                                    thrown = thrown.InnerException;
                                }
                            }
                            reads.CompleteAdding();
                        });                  


                    // Consume them into output files.
                    foreach(var r in reads.GetConsumingEnumerable()) {
                        foreach(var outputter in outputters) {
                            outputter.ConsumeCCSRead(r.Item1, r.Item2, r.Item3, treatment, cov_level);
                        }
                    }

                    // throw any exceptions (this should be used after putting the consumer on a separate thread)
                    producer.Wait();
                    // Close the files
                    outputters.ForEach(z => z.Finish());

            }
            }
            catch(DllNotFoundException thrown) {
                Console.WriteLine ("Error thrown when attempting to generate the CCS results.");
                Console.WriteLine("A shared library was not found.  To solve this, please add the folder" +
                    " with the downloaded files libbwasharp and libMonoPosixHelper" +
                    "to your environmental variables (LD_LIBRARY_PATH on Ubuntu, DYLD_LIBRARY_PATH on Mac OS X)."); 
                Console.WriteLine ("Error: " + thrown.Message);
                Console.WriteLine (thrown.StackTrace);

            }
            catch(Exception thrown) {
                Console.WriteLine ("Error thrown when attempting to generate the CCS results");
                Console.WriteLine ("Error: " + thrown.Message);
                while (thrown.InnerException != null) {
                    Console.WriteLine ("Inner Exception: " + thrown.InnerException.Message);
                    thrown = thrown.InnerException;
                }
            }
        }
       static void DisplayHelp() {
            Console.WriteLine ("ccscheck INPUT OUTDIR REF[optional]");
            Console.WriteLine ("INPUT - the input fastq file");
            Console.WriteLine ("OUTDIR - directory to output results into");
            Console.WriteLine ("REF - A fasta file with the references (optional)");
        }                   
    }
}
