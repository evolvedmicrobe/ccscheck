using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;

using Bio;
using Bio.IO;
using Bio.IO.PacBio;
using Bio.IO.FastQ;
using Bio.Core.Extensions;
using Bio.IO.SAM;

namespace ccscheck
{
    public class FastQCCSReader : ISequenceParser
    {

        public IEnumerable<ISequence> Parse (Stream stream)
        {
            FastQParser fqp = new FastQParser ();
            foreach (var seq in fqp.Parse (stream)) {
                var name = seq.ID;
                var sp = name.Split ('/');
                var movie = sp [0];
                var hole = sp [1];
                SAMAlignedSequence sam = new SAMAlignedSequence ();
                sam.QuerySequence = seq;
                sam.OptionalFields.Add (new SAMOptionalField () { Tag = "sn", Value = "f,0,0,0,0" });
                sam.OptionalFields.Add (new SAMOptionalField () { Tag = "rs", Value = "f,0,0,0,0,0,0" });
                sam.OptionalFields.Add (new SAMOptionalField () { Tag = "zs", Value = "f,0,0,0,0,0,0" });
                PacBioCCSRead read = new PacBioCCSRead (sam) {
                    AvgZscore = Single.NaN,
                    HoleNumber = Convert.ToInt32 (hole),
                    Movie = movie
                };
                yield return read;
            }
        }
        public FastQCCSReader ()
        {
        }

        public IAlphabet Alphabet {
            get {
                throw new NotImplementedException ();
            }

            set {
                throw new NotImplementedException ();
            }
        }

        public string Description {
            get {
                throw new NotImplementedException ();
            }
        }

        public string Name {
            get {
                throw new NotImplementedException ();
            }
        }

        public string SupportedFileTypes {
            get {
                throw new NotImplementedException ();
            }
        }



        public ISequence ParseOne (Stream stream)
        {
            throw new NotImplementedException ();
        }
    }
}


