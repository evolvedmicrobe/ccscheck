#include <string>
namespace PacBio {
namespace Variant {



class VariantCaller {
	public:
		VariantCaller(std::string fastaFile, std::string outputFile);
		~VariantCaller();
		static void CallCCSVariants( std::string str);
};

VariantCaller CreateVariantCaller(std::string fastaFile, std::string outputFile);


} // namespace Consensus
} // namespace PacBio