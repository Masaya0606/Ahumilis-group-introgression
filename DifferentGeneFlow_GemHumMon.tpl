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
3
//Migration matrix 0
0 MIG_01 MIG_02
MIG_10 0 MIG_12
MIG_20 MIG_21 0
//Migration matrix 1
0 MIG2_01 MIG2_02
MIG2_10 0 MIG2_12
MIG2_20 MIG2_21 0
//Migration matrix 2
0 0 0
0 0 0
0 0 0
//historical event: time, source, sink, migrants, new deme size, growth rate, migr mat index
8 historical event
CHANGM 0 0 1 RSANC1_gem 0 0
CHANGM 1 1 1 RSANC1_hum 0 0
CHANGM 2 2 1 RSANC1_mon 0 0
TDIV1 0 0 1 RSANC2_gem 0 1
TDIV1 1 1 1 RSANC2_hum 0 1
TDIV1 2 2 1 RSANC2_mon 0 1
TDIV2 0 1 1 RESIZE0$ 0 2
TDIV3 1 2 1 RESIZE1$ 0 2
//Number of independent loci [chromosome] 
1 0
//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
1
//per Block:data type, number of loci, per gen recomb and mut rates
FREQ 1 0.4773642 2e-8
