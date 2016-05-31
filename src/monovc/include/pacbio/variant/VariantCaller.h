#pragma once

#include <string>
#include <pacbio/consensus/Integrator.h>
#include <pacbio/consensus/Polish.h>

namespace PacBio {
namespace Variant {



class VariantCaller {
	public:
		VariantCaller(std::string fastaFile, std::string outputFile);
		~VariantCaller();
		static void CallCCSVariants(PacBio::Consensus::MonoMolecularIntegrator* ai, std::string movieName, long zmw);
};

VariantCaller CreateVariantCaller(std::string fastaFile, std::string outputFile);


} // namespace Consensus
} // namespace PacBio