import os
import pandas as pd


#swarm sample Bedpe and Single Bedpe Overlap
def generate_bedpe_bedpe_overlap_swarm(sample_list, sample_nh_dir, bedpe, otype, overlap_dir, script_dir):

    bedpe_name = os.path.splitext(os.path.basename(bedpe))[0]

    with open(f"{script_dir}/all_{bedpe_name}_overlap.swarm",'w') as file_handler:
        for sample in sample_list:
            file_handler.write(f"bedtools pairtopair -a {sample_nh_dir}/{sample}.bedpe -b {bedpe} -type {otype} > \
{overlap_dir}/{sample}_{bedpe_name}_overlap.txt\n")

    file_handler.close()
    os.system(f"swarm -f {script_dir}/all_{bedpe_name}_overlap.swarm --module bedtools --g 50")
    
#single Bedpe and Bed Overlap
def run_bedpe_bed_overlap(sample, sample_nh_dir, bed, otype, overlap_dir):
    bed_name = os.path.splitext(os.path.basename(bed))[0]
    os.system(f"module load bedtools; \
              bedtools pairtobed -a {sample_nh_dir}/{sample}.bedpe -b {bed} -type {otype} > {overlap_dir}/{sample}_{bed_name}_overlap.txt")
    
#copy the loop files and rename
def get_sample_loop_files(juicer_dir, sample_dir, sample_list):
    for sample in sample_list:
        os.system(f"cp {juicer_dir}/{sample}/aligned/inter_30_loops/merged_loops.bedpe {sample_dir}")
        os.system(f"mv {sample_dir}/merged_loops.bedpe {sample_dir}/{sample}.bedpe")

#reformat the loop files and relocate them
def get_no_header_loop_files(sample_dir, sample_nh_dir, sample_list):
    for sample in sample_list:
        #read the bedpe and prepend 'chr' to the two chromosome cols
        sample_bedpe = pd.read_csv(f"{sample_dir}/{sample}.bedpe",sep='\t')
        sample_bedpe['chr1'] = 'chr' + sample_bedpe['chr1'].astype(str)
        sample_bedpe['chr2'] = 'chr' + sample_bedpe['chr2'].astype(str)
        sample_bedpe=sample_bedpe[['chr1','x1','x2','chr2','y1','y2']]

        sample_bedpe.to_csv(f"{sample_nh_dir}/{sample}.bedpe", sep='\t', header=False, index=None, mode='w+')

#get a df of overlap data between a bedpe and bed 
def get_bedpe_bed_overlap_data(overlap_file, sample_file):
    
    sample_name = os.path.splitext(os.path.basename(sample_file))[0]
    sample_lines = len(pd.read_csv(sample_file,sep="\t",header=None).index)
    #read the overlap output
    #overlap_data = pd.read_csv(f"{overlap_dir}/{sample}_{bed_name}_overlap.txt",sep="\t",header=None)
    if(os.path.getsize(overlap_file) == 0):
    	return pd.DataFrame(data={'sample':[sample_name], 'counts':0, 'total':[sample_lines], 'percent':[0]})
    overlap_data = pd.read_csv(overlap_file, sep = "\t", header = None)
    
    #we only want the first 6 cols corresponding to the bedpe loops for now
    overlap_subset = overlap_data.iloc[:,:6]
    #remove duplicate rows
    unique_overlaps = overlap_subset.drop_duplicates()
    #count the number of overlaps
    overlap_counts = len(unique_overlaps.index)
    
    
    
    ret_df = pd.DataFrame(data={'sample':[sample_name], 'counts':[overlap_counts], 'total':[sample_lines], 'percent':[(overlap_counts/sample_lines * 100)]})
    
    return ret_df

#All sample bedpe overlap with one bed file
def generate_bedpe_bed_overlap_swarm(sample_list, sample_nh_dir , bed, otype, overlap_dir, script_dir):

    bed_name = os.path.splitext(os.path.basename(bed))[0]

    with open(f"{script_dir}/all_{bed_name}_overlap.swarm",'w') as file_handler:
        for sample in sample_list:
            file_handler.write(f"bedtools pairtobed -a {sample_nh_dir}/{sample}.bedpe -b {bed} -type {otype} > \
{overlap_dir}/{sample}_{bed_name}_overlap.txt\n")

    file_handler.close()
    !swarm -f {script_dir}/all_{bed_name}_overlap.swarm --module bedtools --g 50
    
    
    
#get a df of overlap data between a bedpe and bedpe
def get_bedpe_bedpe_overlap_data(overlap_file, sample_file):
    
    sample_name = os.path.splitext(os.path.basename(sample_file))[0]
    sample_lines = len(pd.read_csv(sample_file,sep="\t",header=None).index)
    #read the overlap output
    #overlap_data = pd.read_csv(f"{overlap_dir}/{sample}_{bed_name}_overlap.txt",sep="\t",header=None)
    if(os.path.getsize(overlap_file) == 0):
    	return pd.DataFrame(data={'sample':[sample_name], 'counts':0, 'total':[sample_lines], 'percent':[0]})
    overlap_data = pd.read_csv(overlap_file, sep = "\t", header = None)
    
    #we only want the first 6 cols corresponding to the bedpe loops for now
    overlap_subset = overlap_data.iloc[:,:6]
    #remove duplicate rows
    unique_overlaps = overlap_subset.drop_duplicates()
    #count the number of overlaps
    overlap_counts = len(unique_overlaps.index)
    
    
    
    ret_df = pd.DataFrame(data={'sample':[sample_name], 'counts':[overlap_counts], 'total':[sample_lines], 'percent':[(overlap_counts/sample_lines * 100)]})
    
    return ret_df

#get a df for all sample bedpe overlap with another bedpe
def get_bedpe_list_bedpe_overlap_data(sample_list, sample_nh_dir , bedpe, overlap_dir, delete = False):

    overlap_df = pd.DataFrame()
    bedpe_name = os.path.splitext(os.path.basename(bedpe))[0]
    for sample in sample_list:
        overlap_file = overlap_dir + '/'+sample+'_'+bedpe_name+'_overlap.txt'
        sample_file = sample_nh_dir + '/' + sample + '.bedpe'
        sample_overlap = get_bedpe_bedpe_overlap_data(overlap_file = overlap_file, sample_file = sample_file)
        overlap_df = overlap_df.append(sample_overlap, ignore_index = True)
        print(sample)
        if(delete):
            os.remove(overlap_file)
    return overlap_df

#get a df for all sample bedpe overlap with a bed file
def get_bedpe_list_bed_overlap_data(sample_list, sample_nh_dir , bed, overlap_dir, delete = False):

    overlap_df = pd.DataFrame()
    bed_name = os.path.splitext(os.path.basename(bed))[0]
    for sample in sample_list:
        overlap_file = overlap_dir + '/'+sample+'_'+bed_name+'_overlap.txt'
        sample_file = sample_nh_dir + '/' + sample + '.bedpe'
        sample_overlap = get_bedpe_bed_overlap_data(overlap_file = overlap_file, sample_file = sample_file)
        overlap_df = overlap_df.append(sample_overlap, ignore_index = True)
        print(sample)
        if(delete):
            os.remove(overlap_file)
    return overlap_df


#overlap samples with each other
def generate_bedpe_between_overlap_swarm(sample_list, sample_nh_dir, otype, overlap_dir, script_dir):

    with open(f"{script_dir}/all_samples_between_overlap.swarm",'w') as file_handler:
        for sample1 in sample_list:
        
            for sample2 in sample_list:
                if(sample1 != sample2):
                    file_handler.write(f"bedtools pairtopair -a {sample_nh_dir}/{sample1}.bedpe -b {sample_nh_dir}/{sample2}.bedpe -type {otype} > \
{overlap_dir}/{sample1}_{sample2}_overlap.txt\n")
            


    file_handler.close()
    os.system(f"swarm -f {script_dir}/all_samples_between_overlap.swarm --module bedtools --g 50")

#get df/matrix of between sample overlap data
def get_bedpe_between_overlap_data(sample_list, sample_nh_dir, overlap_dir, delete = False):

    sample_list = sorted(sample_list)
    all_percent_data={}
    all_count_data={}
    for sample in sample_list:
        if(os.stat(f"{sample_nh_dir}/{sample}.bedpe").st_size != 0):
            sample_lines = len(pd.read_csv(f"{sample_nh_dir}/{sample}.bedpe",sep="\t").index)
        else:
            sample_lines=0

        comp_percent_data={}
        comp_count_data={}

        for comp_sample in sample_list:
            if sample!=comp_sample:
                if(os.stat(f"{overlap_dir}/{sample}_{comp_sample}_overlap.txt").st_size !=0 ):
                    
                    sample_file = f"{sample_nh_dir}/{sample}.bedpe"
                    overlap_file = f"{overlap_dir}/{sample}_{comp_sample}_overlap.txt"
                    df = get_bedpe_bedpe_overlap_data(overlap_file,sample_file)
                    comp_count_data[comp_sample] = df.iloc[0]['counts']
                    comp_percent_data[comp_sample] = df.iloc[0]['percent']
                    


            if sample==comp_sample:
                comp_percent_data[comp_sample]=None
                comp_count_data[comp_sample]=None
        all_percent_data[sample]=comp_percent_data
        all_count_data[sample]=comp_count_data

    percent_df = pd.DataFrame(data=all_percent_data)
    count_df = pd.DataFrame(data=all_count_data)
    return count_df, percent_df

#shuffle a bedpe file n times
def shuffle_bedpe(sample, n, sample_nh_dir, shuffle_dir, sizes):

    #just shuffle once
    if(n==1):
        !(module load bedtools; \
     bedtools shuffle -i {sample_nh_dir}/{sample}.bedpe -g {sizes} -bedpe > {shuffle_dir}/{sample}_shuffle.bedpe)
    else:
        for i in range(n):
            !(module load bedtools; \
     bedtools shuffle -i {sample_nh_dir}/{sample}.bedpe -g {sizes} -bedpe > {shuffle_dir}/{sample}_shuffle_{i}.bedpe)

#shuffle all the sample bedpes n times
def shuffle_bedpe_list(sample_list, n, sample_nh_dir, shuffle_dir, sizes):

    for sample in sample_list:
        shuffle_bedpe(sample = sample, n = n, sample_nh_dir = sample_nh_dir, shuffle_dir = shuffle_dir, sizes = sizes)

#check for overlap with a bed file in all shuffled sample bedpes
def shuffle_swarm_overlap(sample_list, n, sample_nh_dir , bed, otype, shuffle_dir, overlap_dir, script_dir):

    bed_name = os.path.splitext(os.path.basename(bed))[0]
    

    with open(f"{script_dir}/all_{bed_name}_shuffle_overlap.swarm",'w') as file_handler:
        if(n==1):
            for sample in sample_list:
                file_handler.write(f"bedtools pairtobed -a {shuffle_dir}/{sample}_shuffle.bedpe -b {bed} -type {otype} > \
{overlap_dir}/{sample}_shuffle_{bed_name}_overlap.txt\n")
        if(n>1):
            for sample in sample_list:
                for i in range(n):
                    file_handler.write(f"bedtools pairtobed -a {shuffle_dir}/{sample}_shuffle_{i}.bedpe -b {bed} -type {otype} > \
{overlap_dir}/{sample}_shuffle_{i}_{bed_name}_overlap.txt\n")
    file_handler.close()
    !swarm -f {script_dir}/all_{bed_name}_shuffle_overlap.swarm --module bedtools --g 50
    
#get overlap data for all shuffled bedpe files
def get_bedpe_list_bed_shuffle_overlap_data(sample_list, shuffle_dir, n, bed, overlap_dir, delete = False):
    overlap_df = pd.DataFrame()
    bed_name = os.path.splitext(os.path.basename(bed))[0]
 
    if(n==1):
        for sample in sample_list:
            overlap_file = overlap_dir + '/'+sample+'_shuffle_'+bed_name+'_overlap.txt'
            sample_file = shuffle_dir + '/' + sample + '_shuffle.bedpe'
            sample_overlap = get_bedpe_bed_overlap_data(overlap_file = overlap_file, sample_file = sample_file)
            overlap_df = overlap_df.append(sample_overlap, ignore_index = True)
            print(sample)
            if(delete):
                os.remove(overlap_file)
                os.remove(sample_file)
    elif (n>1):
        for sample in sample_list:
            for i in range(n):
                overlap_file = overlap_dir + '/'+sample+'_shuffle_' + str(i) +'_'+bed_name+'_overlap.txt'
                sample_file = shuffle_dir + '/' + sample + '_shuffle_' + str(i) + '.bedpe'
                sample_overlap = get_bedpe_bed_overlap_data(overlap_file = overlap_file, sample_file = sample_file)
                overlap_df = overlap_df.append(sample_overlap, ignore_index = True)
                print(sample + ' ' + str(i))
                if(delete):
                    os.remove(overlap_file)
                    os.remove(sample_file)
    return overlap_df
