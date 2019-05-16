mkdir results_example

name_inst="karate"
instance="real_networks/"$name_inst".paj"
res_base="results_example/part_"$name_inst
time_base="results_example/time_"$name_inst
pareto_file="results_example/pareto_"$name_inst".txt"

MO=11
IT=50
p=0.1
NF=11
tau=0.5
NG=50
NP=5
NO=0.4

./spectral_clustering.out $MO $IT $p $instance $res_base $time_base $pareto_file $NF $tau $NG $NP $NO

