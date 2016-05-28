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

// Fields for items we want to get
MonoClassField POS;
MonoClassField Ref;
MonoClassField CIGAR;
MonoClassField RevComplement;

void create_mono_runtime() {
	char toConvert[] = "embededccscheck.exe";
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
	mAligner = mono_class_from_name( mImage, "embededccscheck", "Aligner");
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
	char* p;
	MonoObject* exception;
	const char* local_str = str.c_str();
	args[0] = mono_string_new(mDomain, local_str);
	MonoObject* result = mono_runtime_invoke(mAlign, NULL, args, &exception);
	if(exception){
		mono_print_unhandled_exception (exception);
		error("C# Based Error");
	}
	MonoString* aln = (MonoString*) result;
	p = mono_string_to_utf8 (aln);
	std::string toprint(p);
	mono_free(p);
	std::cout << toprint << std::endl;
}

int main() {
	create_mono_runtime();
	VariantCaller::CreateVariantCaller("/Users/nigel/BroadBundle/human_g1k_v37.fasta");
	VariantCaller::CallVariants("gatgggaccttgtggaagaagaggtgccaggaatatgtctgggaagggga");

}