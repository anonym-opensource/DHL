# Fast Maintenance of 2-hop Labels for Shortest Distance Queries on Fully Dynamic Graphs

The introduction of these files are as follows. 


# Data

The data files are zipped together, and should be unzipped first. There are 6 datasets in total, and there are 3 files for each dataset. For example, the Amazon dataset contains the following files: 
1. amazon_graph.txt: the readable contents for the amazon graph, with Jaccard and random edge weights inside.
2. amazon_Jaccard.bin: the binary file of the amazon graph with Jaccard edge weights (the reason for generating binary files is that it is much faster to read binary files than to read raw data files).
3. amazon_random.bin: the binary file of the amazon graph with random edge weights.



# C++ codes 

The C++ source codes for the experiments are in <b>exp.cpp</b>. 

Running these codes requires including some header files in cppheaders_202*****.zip and the Boost library (https://www.boost.org/). 

After making the header files and data files ready, <b>we can first run "generate_L_PPR()" in "exp.cpp" to generate all the initial shortest distances indexes, and then run "exp()" in "exp.cpp" to conduct all the experiments in the paper.</b> Make sure there is enough memory (1 TB RAM). 

Specifically, in the terminal on a Linux server, we can compile and run the above codes using the following commands:
```
g++ -std=c++17 -I/home/boost_1_75_0 -I/root/cppheaders_202***** exp.cpp -lpthread -O3 -o A.out
./A.out
```
, where "-I/home/boost_1_75_0" is to add the path of the boost folder when compiling, "-I/root/cppheaders_202*****" is to add the path of the cppheader folder when compiling, "-lpthread" is for parallel computation, and "-O3" is for compiler optimisation.

To read these C++ codes in detail, it is recommended to start from "exp()", and then go to "exp_element1()". More detailed codes in other regions can then be traced. In particular,
- "build_in_progress/HL/dynamic/PLL_dynamic.h" contains codes of PLL to generate L and PPR.
2. "build_in_progress/HL/dynamic/WeightIncreaseMaintenance_improv_multiThread.h" contains codes of the proposed weight increase maintenance algorithm (FastInM).
3. "build_in_progress/HL/dynamic/WeightDecreaseMaintenance_improv_multiThread.h" contains codes of the proposed weight decrease maintenance algorithm (FastDem).
4. "build_in_progress/HL/dynamic/WeightIncrease2021_multiThread.h" contains codes of the existing InAsyn algorithm.
5. "build_in_progress/HL/dynamic/WeightDecrease2021_multiThread.h" contains codes of the existing DeAsyn algorithm.
6. "build_in_progress/HL/dynamic/WeightDecrease2014_multiThread.h" contains codes of the existing DePLL algorithm.
7. "build_in_progress/HL/dynamic/WeightIncrease2019_multiThread.h" contains codes of the existing InPLL algorithm.
