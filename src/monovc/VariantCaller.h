#include <string>
class VariantCaller {
	public:
		static void CreateVariantCaller(std::string fastaFile);
		static void DestroyVariantCaller();
		static void CallVariants(std::string str);
};