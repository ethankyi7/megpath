//NMFMatrix.h
//Matthew Dyer
//Created on 6/1/2017
//Last Modified: 6/12/2017

#ifndef NMFMATRIX__H
#define NMFMATRIX__H

#include <iostream>
#include <fstream>
#include <Eigen/Core>
#include "ProbFunc.h"
#include "FixedPF.h"
#include "PointerSupport.h"
#include "Range.h"

using namespace std;
using namespace Eigen;

typedef Matrix<ProbFunc*, Dynamic, Dynamic> MatrixXp;

class NMFMatrix{
	public:
		NMFMatrix();
		~NMFMatrix();
		NMFMatrix(int rowss, int cols,ProbFunc* p);
		void read(double* data);
		int read(double* data, Range r);
		void write(double* data);
		int write(double* data, Range r);
		void transition();
		int average(double* data, Range r, double alpha);
		void reset();
		void write(string filename);
		void calculateSize();
		int size(){return bufferedSize;};
		int size(Range r);
		void resize(int newRows, int newCols);
		void shrink(Range r);
		void fixRange(Range r);
		void observeRange(Range r);

		bool isConstrained(int row);
		void createBuffers();
		
		int bufferedSize;
		double* sendBuffer;
		double* recvBuffer;
		int rows(){return matrix.rows();};
		int columns(){return matrix.cols();};
		MatrixXd matrix;
		MatrixXp functions;

		ProbFunc* prototype;
};

#endif
