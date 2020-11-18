I gratefully acknowledge that parts of this program were written by Lazaros Gallos in an implementation of the "MEMB" algorithm and were obtained from: https://hmakse.ccny.cuny.edu/software-and-data/

The program has been edited and augmented by Sarah Thomson to compute "probabilistic" multifractal dimension spectra for a local optima network (based on network edge distances and on edge weights/probabilities). 

The algorithm is an updated variant based upon the general template of algorithm 1 in [1] and will be described in the upcoming publication "The Fractal Geometry of Fitness Landscapes at the Local Optima Level." It is a specialised implementation of the "Sandbox" algorithm, which was proposed for fractal analysis of complex networks in [2].

To compile: gcc pf1.c -o multifractal-analysis -lm

... or "make" with the makefile provided.

To run you will need a local optima network in Pajek format including fitness values for each node and standardised edge weights. Nodes must be named as 0 - (n-1). Each line of the input network text file will be: NODENAME1 NODEFITNESS1 NODENAME2 NODEFITNESS2 EDGEWEIGHT\n 

To run: ./multifractal-analysis INPUTNETWORK.txt N OUTPUTFILE.txt NETWORKDIAMETER NUMBERCENTRES DISTANCESTABLE.TXT WEIGHTEDADJACENCIES.TXT

where N is the number of network nodes; NETWORKDIAMETER is the diameter; NUMBERCENTRES is the number of sandbox centres, DISTANCESTABLE.TXT is a text file containing a matrix which is N*N (N being network size) containing the pairwise distances between nodes; WEIGHTEDADJACENCIES.TXT is a text file containing an N*N matrix with pairwise edge weights (if there is an edge) between nodes.

An example: ./multifractal-analysis had12-l0-p3.txt 133 had12-l0-p3fractal.txt 6 50 had12-l0-p3.distancetable had12-l0-p3.weightedadjacency

[1] Thomson, S. L., Verel, S., Ochoa, G., Veerapen, N., & Cairns, D. (2018, July). Multifractality and dimensional determinism in local optima networks. In Proceedings of the Genetic and Evolutionary Computation Conference (pp. 371-378).

[2] Liu, Jin-Long, Zu-Guo Yu, and Vo Anh. "Determination of multifractal dimensions of complex networks by means of the sandbox algorithm." Chaos: An Interdisciplinary Journal of Nonlinear Science 25.2 (2015): 023103.
