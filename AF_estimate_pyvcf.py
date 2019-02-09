#!/usr/bin/env python
# -*- coding: ASCII -*-

#####This python script takes a vcf file and calculates allele frequencies using genotype likelihoods for a subset of individuals at specific loci using the method of Kim et al. BMC Bioinformatics 2011, 12:231 Estimation of allele frequency and associationmapping using next-generation sequencing data
#####Only biallelic loci with a reference allele and one alternate allele.are considered. 
#####numpy, scipy and pyvcf must be installed in the version of python used.
#####vcf files must be bgzipped and tabix indexed.
#####The sample_pop_file file must have the tab separated fields in the following order: sample_name\tpopulation.
#####Each sample should be on a different line, samples can appear in more than one population.
#####The snp_position_file must have the tab separated fields in the following order: chromosome\tsnp_position.
#####Each snp should be on a different line.
#####95% confidence intervals are calculated the using likelihood profile method.
#####Usage is ./AF_estimate_pyvcf.py <in.vcf.gz> <sample_pop_file> <snp_position_file> <min_number_of_samples_per_population_to_consider>
#####Output will be to screen and also written to the sample_pop_file name appended with '.AF;
#####Written (poorly) by Krishna Veeramah (krishna.veeramah@stonybrook.edu)

import string
import os
from sys import argv
import numpy as np
import gzip
from scipy.stats import binom
from scipy.misc import comb
import scipy
import vcf

invcf=argv[1] #tabix indexed vcf
infile=argv[2] #list of samples and populations, tab seperated
snp_pos=argv[3] #list of snps to estimate allele frequencies
min_val=int(argv[4]) #minimum sample size to report

precision=3
e_use=10**-precision

file=open(infile,'r')
data=file.read()
data=string.split(data,'\n')
if data[-1]=='':
    del(data[-1])

pops={}
pops_info={}
for g in range(len(data)):
    k=string.split(data[g],'\t')
    if pops.has_key(k[1])==False:
        pops[k[1]]=[]
        pops_info[k[1]]=k[1:]
    pops[k[1]].append(k[0])


file=open(snp_pos,'r')
data=file.read()
data=string.split(data,'\n')
if data[-1]=='':
    del(data[-1])

snps=[]
for g in range(len(data)):
    k=string.split(data[g],'\t')
    chrom=k[0]
    end=int(k[1])
    start=end-1
    snps.append([chrom,start,end])


#likelihood estimator of allele frequency estimation function from Kim et al. BMC Bioinformatics 2011, 12:231 Estimation of allele frequency and associationmapping using next-generation sequencing data
def AF_estimate(p,GL_array):  #p is frequency of major allele
    p_array=np.array([p**2,2*p*(1-p),(1-p)**2])
    ll=np.prod(np.sum(GL_array*p_array,axis=1))
    return -ll

#1st derivate of allele frequency estimation function
def AF_estimate_der(p,GL_array):  
    p_array=np.array([p**2,2*p*(1-p),(1-p)**2])
    ind_ll=np.sum(GL_array*p_array,axis=1)
    p_array_der=np.array([2*p,2-4*p,2*(p-1)])
    indloci_der=np.sum(GL_array*p_array_der,axis=1)
    der1=np.sum((indloci_der*(np.prod(ind_ll)))/ind_ll)
    return -der1

#likelihood estimator of allele frequency estimation function for scipy
def AF_estimate_ll(p,GL_array):  
    p_array=np.array([p[0]**2,2*p[0]*(1-p[0]),(1-p[0])**2])
    ll=np.prod(np.sum(GL_array*p_array,axis=1))
    return -ll

#1st derivate of allele frequency estimation function for scipy
def AF_estimate_der_ll(p,GL_array):  
    p_array=np.array([p[0]**2,2*p[0]*(1-p[0]),(1-p[0])**2])
    ind_ll=np.sum(GL_array*p_array,axis=1)
    p_array_der=np.array([2*p[0],2-4*p[0],2*(p[0]-1)])
    indloci_der=np.sum(GL_array*p_array_der,axis=1)
    der1=np.sum((indloci_der*(np.prod(ind_ll)))/ind_ll)
    return np.array(-der1)

#likelihood estimator and 1st derivate of allele frequency estimation function for scipy
def AF_estimate_fun_and_der_ll(p,GL_array):  
    p_array=np.array([p[0]**2,2*p[0]*(1-p[0]),(1-p[0])**2])
    ind_ll=np.sum(GL_array*p_array,axis=1)
    ll=np.prod(ind_ll)
    p_array_der=np.array([2*p[0],2-4*p[0],2*(p[0]-1)])
    indloci_der=np.sum(GL_array*p_array_der,axis=1)
    der1=np.sum((indloci_der*(np.prod(ind_ll)))/ind_ll)
    return ll,np.array(-der1)

#likelihood profile-based CIs
def AF_estimate_CI(ml_p,GL_array,e=e_use,C=3.841):  #e is required precision, C is chi-square critical value
    ml=-AF_estimate(ml_p,GL_array)

    chi_sq=0.0
    p=ml_p
    if p<=0.0:
        CI5=0.0
    else:
        while (chi_sq<C):
            ll=-AF_estimate(p,GL_array)
            chi_sq=-2*np.log(ll/ml)
            p-=e
            if p<0.0:
                break
    
        CI5=p+(e*2)

    chi_sq=0.0
    p=ml_p
    if p>=1.0:
        CI95=1.0
    else:
        while chi_sq<C:
            ll=-AF_estimate(p,GL_array)
            chi_sq=-2*np.log(ll/ml)
            p+=e
            if p>1.0:
                break

        CI95=p-(e*2)

    return np.array([CI5,CI95])
    

vcf_reader = vcf.Reader(open(invcf, 'r'))

##all sampls in vcf file
samplist=vcf_reader.samples

##find out sample position indexes for each population
pop_index={}
for i in pops:
    samps=pops[i]
    pop_index[i]=[]
    for gg in range(len(samps)):
        try:
            pop_index[i].append(samplist.index(samps[gg]))
        except:
            ok=1
    if len(pop_index[i])==0:
        del(pop_index[i])


fileout=open(infile+'.AF','w')
out='chrom\tpos\tref_allele\talt_allele\talt_allele_freq\t95%CI\tmaximum_likelihood\tcalled_genotypes\tallele_depths\tnumber_of_usable_samples\tother_meta_information\n'
fileout.write(out)

#iterate through SNPs
for g in range(len(snps)):
    for record in vcf_reader.fetch(snps[g][0],snps[g][1],snps[g][2]):
        pos=record.POS
        chrom=record.CHROM
        ref=record.REF
        alts=record.ALT

        if len(alts)>1:
            print 'More than one alternate allele, skipping '+snps[g][0]+':'+str(snps[g][2])
        elif alts=='.':
            print 'Unknown alternate allele, skipping '+snps[g][0]+':'+str(snps[g][2])
        else:
            alt=str(alts[0])

            for i in pop_index:
                samp_index=pop_index[i]
            
                PLs=[]
                GTs=''
                ADs=''
                for gg in range(len(samp_index)):
                    if record.samples[samp_index[gg]]['GT']<>'./.':
                        PLs.append(record.samples[samp_index[gg]]['PL'])
                        GTs=GTs+record.samples[samp_index[gg]]['GT']+':'
                        ADs=ADs+string.join(record.samples[samp_index[gg]]['AD'],',')+':'

                GLs=10**(-np.array(PLs)/10.0)
                
                nb_use=len(PLs)

                if nb_use>=min_val:
                    temp_grid=[]
                    for g in range(0,11):                                                                                                
                        res=AF_estimate(g/10.0,GLs)
                        temp_grid.append([res,g/10.0])
                    temp_grid.sort()

                    ###try bounded L-BFGS-B from multiple start sites, with derivative
                    temp1=scipy.optimize.minimize(AF_estimate_ll,0.5,args=(GLs,),method='L-BFGS-B',bounds=((0.0,1.0),),jac=AF_estimate_der_ll)
                    temp2=scipy.optimize.minimize(AF_estimate_ll,0.0,args=(GLs,),method='L-BFGS-B',bounds=((0.0,1.0),),jac=AF_estimate_der_ll)
                    temp3=scipy.optimize.minimize(AF_estimate_ll,1.0,args=(GLs,),method='L-BFGS-B',bounds=((0.0,1.0),),jac=AF_estimate_der_ll)
                    temp4=scipy.optimize.minimize(AF_estimate_ll,temp_grid[0][1],args=(GLs,),method='L-BFGS-B',bounds=((0.0,1.0),),jac=AF_estimate_der_ll)

                    all_ll=[[temp1['fun'],temp1['x']],[temp2['fun'],temp2['x']],[temp3['fun'],temp3['x']],[temp4['fun'],temp4['x']]]
                    all_ll.sort()

                    ml_p=all_ll[0][1][0]
                    CI=AF_estimate_CI(ml_p,GLs)
                                  
                    if nb_use>=min_val:
                        #out=string.join(pops_info[i],'\t')+'\t'+chrom+'\t'+str(pos)+'\t'+str(nb_use)+'\t'+ref+':'+str(round(ml_p,3))+'\t'+alt+':'+str(round(1.0-ml_p,3))+'\t'+str(round(CI[0],3))+'-'+str(round(CI[1],3))+'\tll='+str(all_ll[0][0])+'\t'+GTs[:-1]+'\t'+ADs[:-1]+'\n'
                        out=chrom+'\t'+str(pos)+'\t'+ref+'\t'+alt+'\t'+str(round(1-ml_p,precision))+'\t'+str(round(1-CI[1],precision))+'-'+str(round(1-CI[0],precision))+'\tll='+str(all_ll[0][0])+'\t'+GTs[:-1]+'\t'+ADs[:-1]+'\t'+str(nb_use)+'\t'+string.join(pops_info[i],'\t')+'\n'

                        print out[:-1]
                        fileout.write(out)

fileout.close()
