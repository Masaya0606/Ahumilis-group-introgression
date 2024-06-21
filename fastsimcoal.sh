#!/bin/bash
#$ -S /usr/bin/bash
#$ -cwd
#$ -l s_vmem=50G
#$ -l mem_req=50G

mkdir ConstantGeneFlow
cd ConstantGeneFlow

for i in {1..100}
 do
   mkdir ConstantGeneFlow_GemHumMon_run"$i"
   cp ../ConstantGeneFlow_GemHumMon.est  ../ConstantGeneFlow_GemHumMon.tpl  ../ConstantGeneFlow_GemHumMon_jointMAFpop1_0.obs ../ConstantGeneFlow_GemHumMon_jointMAFpop2_0.obs  ../ConstantGeneFlow_GemHumMon_jointMAFpop2_1.obs  ConstantGeneFlow_GemHumMon_run"$i""/"
   cd ConstantGeneFlow_GemHumMon_run"$i"
   singularity exec /usr/local/biotools/f/fastsimcoal2\:27093--hdfd78af_0 fsc27093 -t ConstantGeneFlow_GemHumMon.tpl  -e ConstantGeneFlow_GemHumMon.est -m -0 -c 10 -n 10000 -L 40 -s0 -M -q --foldedSFS
   cd ..
 done
cd ../
