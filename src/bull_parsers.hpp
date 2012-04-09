#if !defined(BULLPARSERS_HPP)
#define BULLPARSERS_HPP

#include <string>
#include <vector>

#include "basic_bull.hpp"
#include "char_encoding.hpp"

namespace bull {

class Tree;

class GetTreesOpts
{
	public:
		GetTreesOpts()
			:storeBrLensFromFile(true),
			fromTree(0),
			toTree(-1)
			{}
			
		bool storeBrLensFromFile;
		std::string filename;
		int fromTree; //index of first tree to store
		int toTree; //index of last tree to get or -1 for all trees
		
		int mode; //paup gettree modes
};

class CodLikeStartOpts
{
	public:
		CodLikeStartOpts()
			:genetic_code(MITO),
			currbrlen(false),
			gtrParams(N_MUT_PARAMS,.25),
			treeScale(1.0)
			{}

		GenCode genetic_code;
		bool currbrlen;
		std::vector<double> gtrParams;
		DblMatrix aaFreqs;
		DblVector multipliers;
		double treeScale;
};

class SimulateOpts
{
	public:
		SimulateOpts()
			:seed(0),
			nReps(1),
			collapsedTreeFilename(),
			concatenations(1),
			overwrite(false),
			paupBlockFilename(),
			outputFilenames(),
			automatic(false),
			verbose(false)
			{}
		
		long seed;
		unsigned nReps; ///# of sim replicates
		unsigned nSimChars;
		std::string collapsedTreeFilename; ///name of file that will store trees that only have internal branches for which there were mutations
		unsigned concatenations; /// # of times (in each replicate) to repeat the "gene" that the model refers to.
		bool overwrite; /// true if output files should be replaced.
		std::vector<std::string> paupBlockFilename;
		std::vector<std::string> outputFilenames; /// datatype to write output in
		std::vector<EncodingType> outTypes; /// datatype to write output in
		std::vector<std::string> treeNames; /// names that the user used to specify the trees
		std::vector<Tree *> treeAlias; ///aliases to the model trees to use
		bool automatic; // true to generate the datatypes for each of the outputs automagically
		std::string tagname;
		bool verbose;
};

} //namespace bull

#endif
