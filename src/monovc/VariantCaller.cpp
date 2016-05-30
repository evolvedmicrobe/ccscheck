#include <exception>
#include <iostream>
#include <stdint.h>
#include <stdio.h>      /* printf, fopen */
#include <stdlib.h>

#include <mono/jit/jit.h>
#include <mono/metadata/assembly.h>
#include <mono/metadata/threads.h>

#include <pacbio/variant/VariantCaller.h>
#include "MonoEmbedding.h"

void error(std::string msg) {
	std::cout << msg << std::endl;
	exit(-1);
}


// Variables in global namespace for now
MonoDomain* mDomain;
MonoClass* mAligner;
MonoMethod* mAlign;
MonoMethod* mSetReferenceFastaAndOutputFile;
MonoMethod* mCloseFileStream;
MonoAssembly* mAssembly;
MonoImage* mImage;


void create_mono_runtime() {
	char toConvert[] = "EmbeddedCCSCheck.exe";
	char* filename = toConvert;
	mDomain = mono_jit_init ("variantcaller");

	mAssembly = mono_domain_assembly_open (mDomain, filename);
	if (!mAssembly) {
		error("Mono assembly seems to be NULL.");
		return;
	}
	mImage = mono_assembly_get_image (mAssembly);
	if (!mImage) {
		error("Mono image seems to be NULL.");
		return;
	}
	mono_jit_exec (mDomain, mAssembly, 1, &filename);
	mAligner = mono_class_from_name( mImage, "EmbeddedCCSCheck", "Aligner");
	mAlign = mono_class_get_method_from_name(mAligner, "Align", 5);
	mCloseFileStream = mono_class_get_method_from_name(mAligner, "CloseFileStream", 0);	
	mSetReferenceFastaAndOutputFile = mono_class_get_method_from_name(mAligner, "SetReferenceFastaAndOutputFile", 2);
}

namespace PacBio {
namespace Variant {


VariantCaller::VariantCaller(std::string fastaFile, std::string outputFile) {
	if (mDomain != NULL) {
		error("Mono domain already created but tried to reinitialize.");
	}

	::mono_mkbundle_init();

	create_mono_runtime();
	const char* loc_fasta_str = fastaFile.c_str();
	const char* loc_output_str = outputFile.c_str();
	void* args[2];
	MonoObject* exception;
	args[0] = mono_string_new(mDomain, loc_fasta_str);	
	args[1] = mono_string_new(mDomain, loc_output_str);
	mono_runtime_invoke(mSetReferenceFastaAndOutputFile, NULL, args, &exception);
	if(exception){
		mono_print_unhandled_exception (exception);
		error("Error initializing the aligner and output files.");
	}
}


VariantCaller::~VariantCaller() {
	MonoObject* exception;
	mono_runtime_invoke(mCloseFileStream, NULL, NULL, &exception);
		if(exception){
		mono_print_unhandled_exception (exception);
		error("Could not close file handle.");
	}
	mono_jit_cleanup (mDomain);
}

void
VariantCaller::CallCCSVariants(std::string str) {
	// Very unclear about this business
	mono_thread_attach (mDomain);
	void* args[5];
	MonoObject* exception;
	const char* local_str = str.c_str();
	long ZMW = 25L;

	args[0] = mono_string_new(mDomain, local_str); // Sequence
	args[1] = &ZMW; // Pointer to abstract integrator
	args[2] = mono_string_new(mDomain, local_str); // Movie Name
	args[3] = &ZMW; // ZMW
	args[4] = &ZMW; // NumReads = both scorable and not.
	std::cout << "So far so good" << std::endl;
	mono_runtime_invoke(mAlign, NULL, args, &exception);
	if(exception){
		mono_print_unhandled_exception (exception);
		error("Error in managed code base when calling variants.");
	}
}
} // namespace Consensus
} // namespace PacBio


// Functions available for PInvoke by C# code.
extern "C" {


	void ScoreVariant(void* ai, int pos, int type, char* bases, double* outputArray) {
		std::string toMutate(bases);

		std::cout << "PTR" << ((long)ai) << std::endl;
		
		std::cout << "POS " << pos << std::endl;
		std::cout << "Type " << type << std::endl;
		std::cout << "Bases " << bases << std::endl;
		outputArray[0] = 2.0;
		outputArray[1] = 3.0;
		}

	void GetBaseLineLikelihoods(void* ai, double* outputArray) {
		for(int i=0; i < 10; i++) {
			outputArray[i] = (double)i;
		}
	}
}

// Only do this if we are not being called as a library, for testing code
#ifdef PROGRAM

int main() {
	using namespace PacBio::Variant;
	VariantCaller vc("/Users/nigel/pacbio/AllReferences.fna", "tmp.txt");
	vc.CallCCSVariants("TGATGATATTGAACAGGAAGGCTCTCCCGACGTTCCGGGTGACAAGCGTATTGAAGGCTC");

}
#endif