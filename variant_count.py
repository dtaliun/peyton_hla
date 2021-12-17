#!/usr/bin/env python
# coding: utf-8

# In[2]:


import pysam
import argparse
import pybedtools
from pybedtools import BedTool


# In[5]:


argparser = argparse.ArgumentParser(description = 'Script to sample random genomic intervals')
argparser.add_argument('-i', '--in-vcf', metavar = 'file', dest = 'vcf_names', required = True, nargs = '+', help = 'Input VCF file(s)')
argparser.add_argument('-l', '--length', metavar = 'interval length', dest = 'length', required = True, type = int, help = 'Length of genomic intervals to be generated')
argparser.add_argument('-n', '--n-intervals', metavar = 'number of intervals', dest = 'n_intervals', required = True, type = int, help = 'Number of genomic intervals to be generated')
argparser.add_argument('-s', '--seed', metavar = 'random seed', dest = 'seed', required = True, type = int, help = 'Random seed parameter for randomized list')



# In[3]:

# count_variants function to get:
# 1. total count 
# 2. total count of variants with allele frequency >0.01
# 3. total count of variants with allele frequency >0.05
# 4. total count of singleton variants i.e. variants with allele count (AC) equal to 1

def count_variants(vcf_name, chromosome, start, end):
    with pysam.VariantFile(vcf_name) as vcf:
        count = 0
        count_af1 = 0
        count_af5 = 0
        count_ac1 = 0
        for variant in vcf.fetch(chromosome, start, end):
            if 'PASS' not in variant.filter:
                continue
            if not "AC" in variant.info:
                continue
                            
            assert len(variant.info["AF"]) == 1
            af = variant.info["AF"][0]
            if af > 0.01:
                count_af1 += 1
            if af > 0.05:
                count_af5 += 1
            
            assert len(variant.info["AC"]) == 1
            ac = variant.info["AC"][0]
            if ac == 1:
                count_ac1 += 1
                
            count = count + 1
        
        return (count, count_af1, count_af5, count_ac1)


# In[2]:


if __name__ == '__main__':
    args = argparser.parse_args()
    
    # build a dictionary to know which VCF stores which chromosome
    chrom2vcf = {}
    
    for vcf_name in args.vcf_names:
        with pysam.TabixFile(vcf_name) as vcf:
            for chrom in vcf.contigs:
                chrom2vcf[chrom] = vcf_name
    
    hg19 = pybedtools.chromsizes('hg19')
    x = BedTool()
    random_list = x.random(l=args.length, n=args.n_intervals, genome='hg19', seed=args.seed) # check random seed paramters
    
    # counts = []
    
    with open("output_file_name.txt", "wt") as output_file:
        for row in random_list:
            chr_name = row[0]
            pos1 = int(row[1])
            pos2 = int(row[2])
            print(chr_name, pos1, pos2)

            random_vcf_name = chrom2vcf[chr_name]
            count = count_variants(random_vcf_name, chr_name, pos1, pos2)
            # count, count_af1, count_af5, count_ac1 = count_variants(...)
            #print(count)
            #counts.append(count)
            output_file.write("{}\n".format("\t".join(str(x) for x in count)))
        
    
    
    

