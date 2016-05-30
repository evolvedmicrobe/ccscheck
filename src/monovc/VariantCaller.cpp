#include <iostream>
#include <stdint.h>
#include <stdio.h>      /* printf, fopen */
#include <stdlib.h>

#include "VariantCaller.h"
#include <mono/jit/jit.h>
#include <mono/metadata/assembly.h>
#include <mono/metadata/threads.h>
#include <exception>

void error(std::string msg) {
	std::cout << msg << std::endl;
	exit(-1);
}


// Variables in global namespace for now
MonoDomain* mDomain;
MonoClass* mAligner;
MonoMethod* mAlign;
MonoMethod* mSetReferenceFasta;
MonoAssembly* mAssembly;
MonoImage* mImage;


void create_mono_runtime() {
	char toConvert[] = "EmbeddedCCSCheck.exe";
	char* filename = toConvert;
	mDomain = mono_jit_init ("variantcaller");

	mAssembly = mono_domain_assembly_open (mDomain, filename);
	if (!mAssembly) {
		error("Mono assembly seems to be NULL??");
		return;
	}
	mImage = mono_assembly_get_image (mAssembly);
	if (!mImage) {
		error("Mono image seems to be NULL??");
		return;
	}
	mono_jit_exec (mDomain, mAssembly, 1, &filename);
	mAligner = mono_class_from_name( mImage, "EmbeddedCCSCheck", "Aligner");
	mAlign = mono_class_get_method_from_name(mAligner, "Align", 1);
	mSetReferenceFasta = mono_class_get_method_from_name(mAligner, "SetReferenceFasta", 1);
}

void
VariantCaller::CreateVariantCaller(std::string fastaFile) {
	const char* local_str = fastaFile.c_str();
	void* args[1];
	MonoObject* exception;
	args[0] = mono_string_new(mDomain, local_str);
	mono_runtime_invoke(mSetReferenceFasta, NULL, args, &exception);
	if(exception){
		mono_print_unhandled_exception (exception);
		error("C# Based Error");
	}
}

void
VariantCaller::DestroyVariantCaller() {
	mono_jit_cleanup (mDomain);
}

void
VariantCaller::CallVariants(std::string str) {
	// Very unclear about this business
	mono_thread_attach (mDomain);
	void* args[1];
	MonoObject* exception;
	const char* local_str = str.c_str();
	args[0] = mono_string_new(mDomain, local_str);
	mono_runtime_invoke(mAlign, NULL, args, &exception);
	if(exception){
		mono_print_unhandled_exception (exception);
		error("C# Based Error");
	}

}

// Functions to make it possible for C++ code to score
// the variants found.
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

}


int main() {
	create_mono_runtime();
	VariantCaller::CreateVariantCaller("/Users/nigel/BroadBundle/human_g1k_v37.fasta");
	VariantCaller::CallVariants("gatgggaccttgtggaagaagaggtgccaggaatatgtctgggaagggga");

}