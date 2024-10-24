//Parameters for the coalescence simulation program : simcoal.exe
3 samples to simulate :
//Population effective sizes (number of genes)
NPOP0
NPOP1
NPOP2
//Samples sizes and samples age 
30
8
13
//Growth rates: negative growth implies population expansion
0
0
0
//Number of migration matrices : 0 implies no migration between demes
2
//Migration matrix 0
0 MIG_01 MIG_02
MIG_10 0 MIG_12
MIG_20 MIG_21 0
//Migration matrix 1
0 0 0
0 0 0
0 0 0
//historical event: time, source, sink, migrants, new deme size, growth rate, migr mat index
5 historical event 
Tadmix 0 0 0 RSANC_ge 0 0
Tadmix 2 2 0 RSANC_hu 0 0
Tadmix 1 1 0 RSANC_mon 0 0 
TDIV2 0 1 1 RESIZE0$ 0 1
TDIV3 1 2 1 RESIZE1$ 0 1
//Number of independent loci [chromosome] 
1 0
//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
1
//per Block:data type, number of loci, per gen recomb and mut rates
FREQ 1 0.4773642 2e-8
