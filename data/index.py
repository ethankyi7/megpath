#write the argument file for output.csv
###
#   outfile.open("test_arguments.txt");
# 	outfile << "analysis = \"test\"\n";
# 	outfile << "stats = \"notAll\"\n";
# 	outfile << "max_runs = 40000\n";
# 	outfile << "print_runs = 10000\n";
# 	outfile << "test_patterns = [\"\",\"\",\"\"]\n";
# 	outfile << "test_filename = \"../testing/testnmf/test_multiplied.csv\"\n";
# 	outfile << "test_origin = [0,0]\n";
# 	outfile << "test_directory = \"results\"\n";
# 	outfile << "test_columns = [0,1,2,3,4]\n";
###

with open('pool_arguments.txt', 'w') as file:
    file.write("analysis = \"pool\"\n")
    file.write("stats = \"all\"\n")
    file.write("max_runs = 40000\n")
    file.write("print_runs = 10000\n")
    file.write("pool_filename = \"../data/pool_output.csv\"\n")
    file.write("pool_origin = [0,0] \n")
    file.write("pool_patterns = [\"\",\"\",\"\"]\n")
    file.write("pool_directory = \"pool_results\"\n")
    testCols = [i for i in range(20)]
    file.write("pool_columns = [" + ",".join(str(i) for i in testCols) + "]\n")
