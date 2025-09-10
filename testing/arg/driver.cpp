//test args driver
//Matthew Dyer
//Created on 5/23/2017
//Last Modified: 5/23/2017

#include<iostream>
#include "../../shared/Value.h"
#include "../../shared/ArgFile.h"

using namespace std;

int main(){
	ArgFile args;

	cout << "something worked" << endl;

	args.fromString("analysis = default\nmax_runs = 1000\ndebug = false\nstart_error = 0.2\nend_error = 0.001\nstart_prob = 0.67\nend_prob = 0.1\nstats = none\nanneal_cut_off = 1.5\ndefault_filename = mixed.csv\ndefault_patterns = {0,0,0,0,0}\ndefault_origin = {1,1}\ndefault_directory = ../testing/csv/\ndefault_controls = {0,0,0,0,0}\nprint_runs = 1000\ninterrupt_runs = 1000\ndefault_ids = {0,0,0,0,0}");
	args.load("tests.txt");
	args.load("ovwrt.txt");

	cout << "Print out the args: \n";
	string print = args.toString();
	cout << print << "\n";

	cout << "Is time an argument?\n";
	cout << args.isArgument("time") << "\n\n";

	cout << "Get the argument for time: \n";
	Value val = args.getArgument("time");
	string out = val;
	cout << out << "\n\n";

	cout << "Vector of ints: \n";
	if(args.isArgument("vecInts")){
		val = args.getArgument("vecInts");
	}
	cout << "[ ";
	vector<int> ints = val;
	for(int i = 0; i < ints.size(); ++i){
		int currInt = ints[i];
		if(i != ints.size()-1){
			cout << currInt << ",";
		}else{
			cout << currInt;
		}
	}
	cout << " ]\n\n";

	cout << "Vector of doubles: \n";
	if(args.isArgument("vecDoubles")){
		val = args.getArgument("vecDoubles");
	}
	cout << "[ ";
	vector<double> doubles;
	for(int i = 0; i < doubles.size(); ++i){
		double currDouble = doubles[i];
		if(i != doubles.size()-1){
			cout << currDouble << ",";
		}else{
			cout << currDouble;
		}
	}
	cout << " ]\n\n";

	cout << "Vector of strings: \n";
	if(args.isArgument("vecStrings")){
		val = args.getArgument("vecStrings");
	}
	cout << "[ ";
	vector<string> strings;
	for(int i = 0; i < strings.size(); ++i){
		string currString = strings[i];
		if(i != strings.size()-1){
			cout << currString << ",";
		}else{
			cout << currString;
		}
	}
	cout << " ]\n\n";

	return 0;
}
