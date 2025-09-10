//OMP Monte Anneal functions
//Dakota Martin
//Julian Dymacek
//Created on 6/7/2018

#include "OMPMonteAnneal.h"

OMPMonteAnneal::OMPMonteAnneal(State* st,int nt): MonteAnneal(st){
	numThreads = nt;
	random_device rd;
	ProbFunc::generator.seed(rd());
	uniform = new UniformPF();
	callback = NULL;
}

/*Run a monte carlo markov chain*/
void OMPMonteAnneal::monteCarloThread(int xStart, int xEnd,int yStart,int yEnd){
	ErrorFunctionRow efRow(state);
	ErrorFunctionCol efCol(state);
	
	//For each spot take a gamble and record outcome
	for(int i =0; i < state->MAX_RUNS; i++){
		Range r;
		r.colStart = 0;
		r.colEnd = state->coefficients.columns();
		r.rowStart = yStart;
		r.rowEnd = yEnd;
		monteCarloStep(state->coefficients,&efRow,r);
        	if(state->both){
			#pragma omp barrier
				Range s;
				s.colStart = xStart;
				s.colEnd = xEnd;
				r.rowStart = 0;
				s.rowEnd = state->patterns.rows();
        		monteCarloStep(state->patterns,&efCol,s);
        	}
		#pragma omp barrier
		#pragma omp master
		{
        	if(i % state->printRuns == 0 && callback != NULL){
			callback->montePrintCallback(i);
        	}
			if(i % state->interruptRuns == 0){
				if(callback != NULL){
					callback->monteCallback(i);
				}
			}
		}
		#pragma omp barrier
	}
}


double OMPMonteAnneal::monteCarlo(){
    Stopwatch watch;
	watch.start();
	
	vector<Range> ranges = state->splitRanges(numThreads);

	int id;
	#pragma omp parallel private(id) num_threads(numThreads)
	{
		id = omp_get_thread_num();
		if(state->constrained){
            ranges[id].colEnd = (id == 0 ? state->patterns.columns() : 0);
        }  
		monteCarloThread(ranges[id].colStart,ranges[id].colEnd,ranges[id].rowStart,ranges[id].rowEnd);
	}

	ErrorFunctionRow efRow(state);
	if(state->debug){
		cout << "Final Error: " << efRow.error() << endl;
    		cout << "Error Histogram: " << efRow.errorDistribution(10) << endl;
    		cout << "Total time: " << watch.formatTime(watch.stop()) << endl;
	}
	return efRow.error();
}


void OMPMonteAnneal::annealThread(int xStart, int xEnd,int yStart,int yEnd){
    ErrorFunctionRow efRow(state);
    ErrorFunctionCol efCol(state);
	double t = state->calcT();
	double alpha = state->calcAlpha(t);

    //For each spot take a gamble and record outcome
    for(int i =0; i < 2*state->MAX_RUNS; i++){
		Range r;
		r.colStart = 0;
		r.colEnd = state->coefficients.columns();
		r.rowStart = yStart;
		r.rowEnd = yEnd;
		annealStep(state->coefficients,t,&efRow,r);
        if(state->both){
			#pragma omp barrier
			Range s;
			s.colStart = xStart;
			s.colEnd = xEnd;
			s.rowStart = 0;
			s.rowEnd = state->patterns.rows();
			annealStep(state->patterns,t,&efCol,s);
		}
		#pragma omp barrier
		#pragma omp master
		{
			if(i % state->interruptRuns == 0 && callback != NULL){
				callback->annealCallback(i);
			}

			if(i % state->printRuns == 0 && callback != NULL){
				callback->annealPrintCallback(i);
			}
		}
		#pragma omp barrier
		t *= alpha;
    }
}

double OMPMonteAnneal::anneal(){
	Stopwatch watch;
	watch.start();
	
	vector<Range> ranges = state->splitRanges(numThreads);
	int id;
	#pragma omp parallel private(id) num_threads(numThreads)
	{
		id = omp_get_thread_num();
		if(state->constrained){
			ranges[id].colStart = 0;
			ranges[id].colEnd = (id == 0 ? state->patterns.columns() : 0);
		}
        annealThread(ranges[id].colStart,ranges[id].colEnd,ranges[id].rowStart,ranges[id].rowEnd);
	}

	ErrorFunctionRow efRow(state);
	if(state->debug){
		cout << "Final Error: " << efRow.error() << endl;
		cout << "Error Histogram: " << efRow.errorDistribution(10) << endl;
		cout << "Total time: " << watch.formatTime(watch.stop()) << endl;
	}
	return efRow.error();
}

