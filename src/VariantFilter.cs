using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;

using Bio;
using Bio.IO.PacBio;
using Bio.Variant;
using Bio.VCF;


namespace ccscheck
{
    /// <summary>
    /// Filters any variants that (in 0 based coordinates) starts at the same location as 
    /// a variant in a VCF file (with 1 based coordinates).  Only filtering on position to avoid 
    /// issues translating positions.
    /// </summary>
    public class VariantFilter
    {
        AATree<Variant> searchTree;

        /// <summary>
        /// Initializes a new instance of the <see cref="ccscheck.VariantFilter"/> class.
        /// 
        /// </summary>
        /// <param name="filename">Filename. Positions must be 1 based.</param>
        public VariantFilter (string filename)
        {
            var variantFile = new VCFParser (filename);
            var variants = variantFile.ToList ();
            searchTree = new AATree<Variant> (new VariantComparer ());
            foreach (var v in variants) {
                if (v.Indel) {
                    var va = new IndelVariant (v.Start - 1, 0, String.Empty, IndelType.Deletion, 'N', -1, false);
                    va.RefName = v.Chr;
                    searchTree.Add (va);
                } else if (v.SNP) {
                    
                    var va = new SNPVariant (v.Start - 1, (char)v.Alleles.First ().Bases [0], (char)v.Reference.Bases [0]);
                    va.RefName = v.Chr;
                    searchTree.Add (va);
                } else {
                    Console.WriteLine ("Skipped complex variant in VCF file: " + v.ToString ());
                }
            }
        }

        public bool ContainsVariantPosition(Variant v)
        {
            return searchTree.Contains (v);
        }
    }

    public class VariantComparer : IComparer<Variant> {

        public int Compare(Variant a, Variant b)
        {
            var car1 = String.Compare(a.RefName, b.RefName, StringComparison.Ordinal);
            if (car1 == 0) {
                return a.StartPosition.CompareTo (b.StartPosition);
            }   
            return car1;
        }
    }
}

