
I gratefully acknowledge that parts of this program were written by Lazaros Gallos in an implementation of the "MEMB" algorithm
and were obtained from: https://hmakse.ccny.cuny.edu/software-and-data/

The program has been edited and augmented by Sarah Thomson to compute "deterministic" multifractal dimension spectra for a local optima network 
(based on network edge distance and on fitness distances).

The algorithm details are described in [1].

It is a specialised implementation of the "Sandbox" algorithm, which was proposed for fractal analysis of complex networks in [2].
To compile: gcc mf1.c -o multifractal-analysis -lm

... or "make" with the makefile provided.

To run you will need a local optima network in Pajek format including fitness values for each node.

Nodes must be named as 0 - (n-1). 

Each line of the input network text file will be: NODENAME1 NODEFITNESS1 NODENAME2 NODEFITNESS2\n 

To run: ./multifractal-analysis INPUTNETWORK.txt N OUTPUTFILE.txt NETWORKDIAMETER NUMBERCENTRES DISTANCESTABLE.TXT

where N is the number of network nodes; NETWORKDIAMETER is the diameter; NUMBERCENTRES is the number of sandbox centres, and DISTANCESTABLE.TXT is a text file containing a matrix which is N*N (N being network size) containing the pairwise distances between nodes.

An example: ./multifractal-analysis had12-l0-p3.txt 133 had12-l0-p3fractal.txt 6 50 had12-l0-p3.distancetable

[1] Thomson, Sarah L., Gabriela Ochoa, and Sébastien Verel. "The fractal geometry of fitness landscapes at the local optima level." Natural Computing (2020): 1-17.
[2] Liu, Jin-Long, Zu-Guo Yu, and Vo Anh. "Determination of multifractal dimensions of complex networks by means of the sandbox algorithm." Chaos: An Interdisciplinary Journal of Nonlinear Science 25.2 (2015): 023103.

