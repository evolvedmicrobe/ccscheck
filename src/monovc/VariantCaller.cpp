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

/* Functions to expose as internal mono calls. */
MonoArray* InternalGetReadNames(PacBio::Consensus::MonoMolecularIntegrator* ai){
	std::vector<std::string> readNames = ai->ReadNames();
	MonoArray* arr = (MonoArray*)mono_array_new (mDomain,  mono_get_string_class (), readNames.size());
	for (int i=0; i < readNames.size(); i++) {
		MonoString* name = mono_string_new(mDomain, readNames[i].c_str()); 
		mono_array_setref(arr, i, name);
	}
	return arr;
}

MonoArray* InternalGetReadDirections(PacBio::Consensus::MonoMolecularIntegrator* ai){
	std::vector<std::string> direcs = ai->ReadDirections();
	MonoArray* arr = (MonoArray*)mono_array_new (mDomain,  mono_get_string_class (), direcs.size());
	for (int i=0; i < direcs.size(); i++) {
		MonoString* name = mono_string_new(mDomain, direcs[i].c_str()); 
		mono_array_setref(arr, i, name);
	}
	return arr;
}

int InternalGetEvaluatorCount(PacBio::Consensus::MonoMolecularIntegrator* ai) {
	return ai->NumEvaluators();
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
	return 'F';
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
	mono_add_internal_call("EmbeddedCCSCheck.VariantScorer::InternalGetEvaluatorCount", reinterpret_cast<void*>(InternalGetEvaluatorCount));
	mono_add_internal_call("EmbeddedCCSCheck.VariantScorer::InternalGetReadDirections", reinterpret_cast<void*>(InternalGetReadDirections));
	

	mAligner = mono_class_from_name( mImage, "EmbeddedCCSCheck", "VariantOutputter");
	mAlign = mono_class_get_method_from_name(mAligner, "AlignAndCallVariants", 4);
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
	void* args[4];
	MonoObject* exception;
	// TODO: Need to fix this copy operation
	//std::string seq = GetSequenceFromAI(ai);
	std::string seq = ai->operator std::string();
	const char* local_str = seq.c_str();
	const char* movie_str = movieName.c_str();
	long ZMW = zmw;

	args[0] = mono_string_new(mDomain, local_str); // Sequence
	args[1] = &ai; // Pointer to abstract integrator
	args[2] = mono_string_new(mDomain, movie_str); // Movie Name
	args[3] = &ZMW; // ZMW
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
	
	void GetBaseLineLikelihoods(PacBio::Consensus::MonoMolecularIntegrator* ai, double* outputArray) {
		std::vector<double> scores = ai->LLs();
		std::copy(scores.begin(), scores.end(), outputArray);
	}

	bool ScoreVariant(PacBio::Consensus::MonoMolecularIntegrator* ai, int pos, int type, char* bases, double* outputArray) {
		bool worked = true;
		try {
			std::string sbases(bases);
			MutationType t = static_cast<MutationType>(type); 
			// If we have only one edit, we will test and not apply the mutation
			if (sbases.size() == 1) {
				//char base = t != MutationType::DELETION ? bases[0] : '-';
				char base = bases[0];
				Mutation m(t, pos, base);
				std::vector<double> scores = ai->LLs(m);
				std::copy(scores.begin(), scores.end(), outputArray);
			} else {
				// We have to apply a series of mulitple mutations to test the likelihood.
				// Note: The mutations are applied in "REVERSE order"
				std::vector<Mutation> rev_muts(sbases.size(), Mutation(t, 0, '-'));
				std::transform(sbases.rbegin(), sbases.rend(), rev_muts.begin(), 
					[t, pos](unsigned char bp) { return Mutation(t, static_cast<size_t>(pos), bp);});

				// Let's get those scores
				ai->ApplyMutations(&rev_muts);
				GetBaseLineLikelihoods(ai, outputArray);

				// Now to get the inverse mutations
				std::vector<Mutation> inverse_muts(rev_muts.size(), Mutation(t, 0, '-'));
				std::transform(rev_muts.begin(), rev_muts.end(), inverse_muts.begin(), 
					[](const Mutation mut) { return GetInverseMutation(mut, mut.Base);});
				// And apply them to undo the damage
				ai->ApplyMutations(&inverse_muts);
			}
		} catch(...) {
				worked = false;
		}
		return worked;
	}



	void GetBaseLineZScores(PacBio::Consensus::MonoMolecularIntegrator* ai, double* outputArray) {
		std::vector<double> zscores = ai->ZScores();
		std::copy(zscores.begin(), zscores.end(), outputArray);
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