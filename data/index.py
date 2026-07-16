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

with open('clusterFourLT_arguments.txt', 'w') as file:
    file.write("analysis = \"clusterFourLT\"\n")
    file.write("stats = \"all\"\n")
    file.write("max_runs = 40000\n")
    file.write("print_runs = 10000\n")
    file.write("clusterFourLT_filename = \"../data/cluster_four.csv\"\n")
    file.write("clusterFourLT_origin = [0,0] \n")
    file.write("clusterFourLT_patterns = [\"patternFourLT\",\"\",\"\"]\n")
    file.write("clusterFourLT_directory = \"clustersLT\"\n")
    testCols = [i for i in range(4)]
    file.write("clusterFourLT_columns = [" + ",".join(str(i) for i in testCols) + "]\n")
    file.write("patternFourLT = [0.009090,0.493705,0.017601,0.024007]")
