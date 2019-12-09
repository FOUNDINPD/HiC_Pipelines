#script to convert cLoop loop files to BEDPE format (formats listed below)
#pass in .loop file as argument and output file directory and include significant loops (1 to include them, 0 to not)
#ex: loop2bedpe.py 35233.loop /output 1
#script by Frank Grenn


#BEDPE
#chrom1	start1	end1	chrom2	start2	end2	name	score	strand1	strand2	

#loop
#loopId	ES	FDR	binomal_p-value	distance	hypergeometric_p-value	iva	ivb	poisson_p-value	ra	rab	rb	poisson_p-value_corrected	binomial_p-value_corrected	hypergeometric_p-value_corrected	significant

#want iva, ivb

import sys
import pandas as pd


cLoop = pd.read_csv(sys.argv[1], sep='\t')

output_path= sys.argv[2]
sig=sys.argv[3]
#subset for only significant loops
if sig=="1" :
	print("only include significant loops!")
	cLoop = cLoop.loc[cLoop['significant'] == 1.0]



ivas=cLoop['iva'].str.split(":|-",2,expand=True)
ivbs=cLoop['ivb'].str.split(":|-",2,expand=True)

bedpe_data = {"chrom1":ivas[0],"start1":ivas[1],"end1":ivas[2],"chrom2":ivbs[0],"start2":ivbs[1],"end2":ivbs[2]}

bedpe=pd.DataFrame(data=bedpe_data)
bedpe=bedpe[['chrom1','start1','end1','chrom2','start2','end2']]

bedpe.to_csv(output_path+'/loops.bedpe',sep='\t', header=False, index=False)