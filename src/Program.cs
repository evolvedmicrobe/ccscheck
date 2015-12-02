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

                if (args.Length > 4) {
                    Console.WriteLine ("Too many arguments");
                    DisplayHelp ();
                } else if (args.Length < 2) {
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
                    BWAPairwiseAligner bwa = null;
                    bool callVariants = ref_name != null;
                    if(callVariants) {
                        bwa = new BWAPairwiseAligner(ref_name, false); 
                    }

                    // Produce aligned reads with variants called in parallel.
                    var reads = new BlockingCollection<Tuple<IQualitativeSequence, BWAPairwiseAlignment, List<Variant>>>();
                    Task producer = Task.Factory.StartNew(() =>
                        {
                            try 
                            {
                                Parallel.ForEach(fqp.Parse(fastq_name), z => {
                                    try {
                                        BWAPairwiseAlignment aln = null;
                                        List<Variant> variants = null;
                                        if (callVariants) {
                                            aln = bwa.AlignRead(z) as BWAPairwiseAlignment;
                                            if (aln!=null) {
                                                variants = VariantCaller.CallVariants(aln);
                                                variants.ForEach( p => {
                                                    p.StartPosition += aln.AlignedSAMSequence.Pos;
                                                    p.RefName = aln.Reference;
                                                    });

                                            }
                                        }
                                        var res = new Tuple<IQualitativeSequence, BWAPairwiseAlignment, List<Variant>>(z, aln, variants);
                                        reads.Add(res);
                                    }
                                    catch(Exception thrown) {
                                        Console.WriteLine("FASTQ READ FAIL: " + z.ID);
                                        Console.WriteLine(thrown.Message);
                                    } });
                            } catch(Exception thrown) {
                                Console.WriteLine("Could not parse FASTQ file: " + thrown.Message);
                                while(thrown.InnerException!=null) {
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
                Console.WriteLine (thrown.StackTrace);
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
