#include <exception>
#include <iostream>
#include <stdint.h>
#include <stdio.h>      /* printf, fopen */
#include <stdlib.h>
#include <vector>

#include <mono/jit/jit.h>
#include <mono/metadata/assembly.h>
#include <mono/metadata/environment.h>
#include <mono/metadata/threads.h>
#include <pacbio/consensus/Mutation.h>

#include <pacbio/variant/VariantCaller.h>
#include "MonoEmbedding.h"


using PacBio::Consensus::Mutation;
using PacBio::Consensus::MutationType;

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

// Pre-declare methods
char revComp(char bases);
Mutation GetInverseMutation(const Mutation& mut, char mutBase);

MonoArray* InternalGetReadNames(PacBio::Consensus::MonoMolecularIntegrator* ai){
	std::vector<std::string> readNames = ai->ReadNames();
	MonoArray* arr = (MonoArray*)mono_array_new (mDomain,  mono_get_string_class (), readNames.size());
	for (int i=0; i < readNames.size(); i++) {
		MonoString* name = mono_string_new(mDomain, readNames[i].c_str()); 
		mono_array_setref(arr, i, name);
	}
	return arr;
}

// Code to test whether applying and unapplying a mutation leads to the expected results.
MonoString* GetTemplateAfterMutation(PacBio::Consensus::MonoMolecularIntegrator* ai, int pos, int type, char mutBase) {
	MutationType t = static_cast<MutationType>(type);
	char base = t != MutationType::DELETION ? mutBase : '-';
	Mutation m(t, pos, base);
	ai->ApplyMutation(m);
	Mutation inverse = GetInverseMutation(m, mutBase);
	std::string seq = ai->operator std::string();
	const char* local_str = seq.c_str();
	ai->ApplyMutation(inverse);
	return mono_string_new(mDomain, local_str);
}

char revComp(char bases) {
	switch (bases) {
		case 'A':
		case 'a':
			return 'T';
		case 'C':
		case 'c':
			return 'G';
		case 'G':
		case 'g':
			return 'C';
		case 'T':
		case 't':
			return 'A';
		case '-':
			return '-';
	}
	error("it's a crazy base in a crazy world.");
}

Mutation GetInverseMutation(const Mutation& mut, char mutBase) {
	MutationType inverseType;
	char base;
	if (mut.Type == MutationType::DELETION) {
		inverseType = MutationType::INSERTION;
		base = mutBase;
	} else if (mut.Type == MutationType::INSERTION) {
		inverseType = MutationType::DELETION;
		base = '-';
	} else {
		inverseType = MutationType::SUBSTITUTION;
		base = revComp(mutBase);
	}
	return Mutation(inverseType, mut.Start(), base);
}

MonoArray* InternalGetReadDirections(PacBio::Consensus::MonoMolecularIntegrator* ai){
	std::vector<std::string> readNames = ai->ReadNames();
	MonoArray* arr = (MonoArray*)mono_array_new (mDomain,  mono_get_string_class (), readNames.size());
	for (int i=0; i < readNames.size(); i++) {
		MonoString* name = mono_string_new(mDomain, readNames[i].c_str()); 
		mono_array_setref(arr, i, name);
	}
	return arr;
}

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
	mono_add_internal_call("EmbeddedCCSCheck.VariantScorer::InternalGetReadNames", reinterpret_cast<void*>(InternalGetReadNames));
	mono_add_internal_call("EmbeddedCCSCheck.VariantScorer::GetTemplateAfterMutation", reinterpret_cast<void*>(GetTemplateAfterMutation));
	
	mAligner = mono_class_from_name( mImage, "EmbeddedCCSCheck", "VariantOutputter");
	mAlign = mono_class_get_method_from_name(mAligner, "AlignAndCallVariants", 5);
	mCloseFileStream = mono_class_get_method_from_name(mAligner, "CloseFileStream", 0);	
	mSetReferenceFastaAndOutputFile = mono_class_get_method_from_name(mAligner, "SetReferenceFastaAndOutputFile", 2);
	
}

std::string GetSequenceFromAI(PacBio::Consensus::MonoMolecularIntegrator* ai) {
    std::string result;
    result.resize(ai->TemplateLength());

    for (size_t i = 0; i < ai->TemplateLength(); ++i)
        result[i] = ai->operator[](i);

    return result;
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
VariantCaller::CallCCSVariants(PacBio::Consensus::MonoMolecularIntegrator* ai, std::string movieName, long zmw) {
	// Very unclear about this business
	mono_thread_attach (mDomain);
	void* args[5];
	MonoObject* exception;
	// TODO: Need to fix this copy operation
	//std::string seq = GetSequenceFromAI(ai);
	std::string seq = ai->operator std::string();
	const char* local_str = seq.c_str();
	const char* movie_str = movieName.c_str();
	long ZMW = zmw;
	int scorable = static_cast<int>(ai->LLs().size());

	args[0] = mono_string_new(mDomain, local_str); // Sequence
	args[1] = &ai; // Pointer to abstract integrator
	args[2] = mono_string_new(mDomain, movie_str); // Movie Name
	args[3] = &ZMW; // ZMW
	args[4] = &scorable; // NumReads = both scorable and not.
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

	void ScoreVariant(PacBio::Consensus::MonoMolecularIntegrator* ai, int pos, int type, char* bases, double* outputArray) {
		MutationType t = static_cast<MutationType>(type);
		char base = t != MutationType::DELETION ? bases[0] : '-';
		Mutation m(t, pos, base);
		std::vector<double> scores = ai->LLs(m);
		std::copy(scores.begin(), scores.end(), outputArray);
	}

	void GetBaseLineLikelihoods(PacBio::Consensus::MonoMolecularIntegrator* ai, double* outputArray) {
		std::vector<double> scores = ai->LLs();
		std::copy(scores.begin(), scores.end(), outputArray);
	}
}

// Only do this if we are not being called as a library, for testing code
#ifdef PROGRAM

int main() {
	using namespace PacBio::Variant;
	//VariantCaller vc("/Users/nigel/pacbio/AllReferences.fna", "tmp.txt");
	//vc.CallCCSVariants("TGATGATATTGAACAGGAAGGCTCTCCCGACGTTCCGGGTGACAAGCGTATTGAAGGCTC");

}
#endif