import sys 
from subprocess import call

sample=sys.argv[1]
root = sample.replace('.txt', '')

call(['perl', 'analysis.pl', sample])

out_files = ['.het.est', '.hm_config.est', '.hm_config.in', '.hm_config.out', '.hm_het.est', '.hm_het.in', '.hm_het.out', '.mesh.in']

for out in out_files:
	new_file=root+out
	call(['mv', out, new_file])
