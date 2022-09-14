
WD=$(cd ../src; pwd)
echo $WD
export PYTHONPATH="${PYTHONPATH}:${WD}"
python3 ../src/main_crossovers.py -i ref_genome_alphaimputeout.haplotypes -p example.fam -m example.map -o testexample