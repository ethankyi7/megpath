//State.h
//Julian Dymacek
//Matthew Dyer
//Dakota Martin
//Created: 6/6/2017
//Modified: 6/10/2018

#ifndef STATE__H
#define STATE__H

#include <cmath>
#include "ArgFile.h"
#include "CSVFile.h"
#include "NMFMatrix.h"
#include  <Eigen/Core>
#include "ShiftPF.h"
#include "WeightedPF.h"
#include "PiecewisePF.h"
#include "PNG.h"

using namespace std;
using namespace Eigen;

class State{
	public:
		State();
		int MAX_RUNS;
		NMFMatrix patterns;
		NMFMatrix coefficients;
		MatrixXd expression;
		string directory;
		string analysis;
		string filename;
		string stats;
		string dist;
		bool both;
		bool debug;
		bool constrained;
		bool img;
		bool gray;
		vector<string> patternNames;
		vector<string> ids;
		int printRuns;
		int interruptRuns;
		int zeroes;
		double annealCutOff;
		double errorAvg;
		double errorCount;
		double start_error;
		double end_error;
		double start_prob;
		double end_prob;
		virtual bool load(string argFileName);
		vector<Range> splitRanges(int by);
		Range getRange(int rank);
		double calcT();
		double calcAlpha(double t);
		void patternMatch(MatrixXd& other);
		void reset();
		void MXdToPNG(MatrixXd mat, int r, int c, bool g, const char* name);
        void MXdToPNG2(MatrixXd mat, int r, int c, int g, const char* name);
		void errorToPNG();
        void errorToPNG2(const char* name);
		void unshufflePC();
		void reshufflePC();

	protected:
		PermutationMatrix<Dynamic> rPerm;
		PermutationMatrix<Dynamic> cPerm;
		PermutationMatrix<Dynamic> zCol;
		void normalize();
		void normalizeMatrix(MatrixXd& mat);
		vector<vector<Value> > pixlToVal(Image* png, bool& gray);
		void sortZero();
};

#endif
