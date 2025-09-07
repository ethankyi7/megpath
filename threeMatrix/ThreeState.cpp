#include "ThreeState.h"
#include <math.h>

ThreeState::ThreeState(){
}
bool ThreeState::load(string argFileName){
	bool rv = State::load(argFileName);
	ArgFile args;
	CSVFile file;
	args.load(argFileName);
	int WEIGHTS = 5;
	int PATTERNS = 0;
	vector<Value> origin;
	
	if(args.isArgument(analysis + "patterns")){
		Value val = args.getArgument(analysis + "patterns");
		vector<Value> patternArgs = val;
		PATTERNS = patternArgs.size();
	}

	if(args.isArgument(analysis + "weights")){
		Value val = args.getArgument(analysis + "wieghts");
		WEIGHTS = val;
	}
	
	if(args.isArgument(analysis + "origin")){
		Value val = args.getArgument(analysis + "origin");
		origin = val;
	}else{
		Value val;
		origin.push_back(val);
		origin.push_back(val);
	}
	vector<vector<Value> > csv = file.readCSV(filename);

	//expression matrix
	int originElement = origin[1];
	int ROWS = csv.size() - originElement;
	coefficients.resize(ROWS,WEIGHTS);
	weights.resize(WEIGHTS,PATTERNS);
	cout << "ROWS:" << ROWS << " WEIGHTS:" << WEIGHTS << " PATTERNS:" << PATTERNS << endl;
	return rv;
}
