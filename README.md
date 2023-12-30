# Fast Maintenance of 2-hop Labels for Shortest Distance Queries on Fully Dynamic Graphs

The introduction of these files are as follows. 


# Data

The data files are zipped together, and should be unzipped first. There are 6 datasets in total, and there are 3 files for each dataset. For example, the Amazon dataset contains the following files: 
1. amazon_graph.txt: the readable contents for the amazon graph, with Jaccard and random edge weights inside.
2. amazon_Jaccard.bin: the binary file of the amazon graph with Jaccard edge weights (the reason for generating binary files is that it is much faster to read binary files than to read raw data files).
3. amazon_random.bin: the binary file of the amazon graph with random edge weights.



# C++ codes 

The codes for conducting the experiments are in <b>exp.cpp</b>. 

The codes for the algorithms are in the header files in cppheaders_202*****.zip.

To read these C++ codes in detail, it is recommended to start from "exp()". More detailed codes in other regions can then be traced. In particular, in the cppheaders_202***** folder,
- "build_in_progress/HL/dynamic/PLL_dynamic.h" contains codes of <b>PLL</b> to generate L and PPR.
- "build_in_progress/HL/dynamic/WeightIncreaseMaintenance_improv.h" contains codes of the proposed weight increase maintenance algorithm (<b>FastInM</b>).
- "build_in_progress/HL/dynamic/WeightDecreaseMaintenance_improv.h" contains codes of the proposed weight decrease maintenance algorithm (<b>FastDeM</b>).
- "build_in_progress/HL/dynamic/WeightIncrease2021.h" contains codes of the existing <b>InAsyn</b> algorithm.
- "build_in_progress/HL/dynamic/WeightDecrease2021.h" contains codes of the existing <b>DeAsyn</b> algorithm.
- "build_in_progress/HL/dynamic/WeightDecrease2014.h" contains codes of the existing <b>DePLL</b> algorithm.
- "build_in_progress/HL/dynamic/WeightIncrease2019.h" contains codes of the existing <b>InPLL</b> algorithm.



# Settings

To compile and run exp.cpp, prepare the environment as follows:

- downlaod and unzip the header files in cppheaders_202*****.zip
- downlaod the Boost library at https://www.boost.org
- downlaod and unzip the datasets
- get a server with enough memory (1 TB RAM) and hard disk space (1 TB). (to only run the small CondMat and Gnutella datasets, just personal computers are OK)

After preparing the environment as suggested above, in the terminal on a Linux server, we can compile and run exp.cpp using the following commands:
```
g++ -std=c++17 -I/home/boost_1_75_0 -I/root/cppheaders_202***** exp.cpp -lpthread -O3 -o A.out
./A.out
```
, where "-I/home/boost_1_75_0" is to add the path of the boost folder when compiling, "-I/root/cppheaders_202*****" is to add the path of the cppheader folder when compiling, "-lpthread" is for parallel computation, and "-O3" is for compiler optimisation.

Specifically, <b>we can first run "generate_L_PPR()" in "exp.cpp" to generate all the initial shortest distances indexes, and then run "exp()" in "exp.cpp" to conduct all the experiments in the paper.</b> 

All the experiments in the paper are conducted on a Linux server with Ubuntu system, 1 TB RAM and 2 TB hard disk space.

