import glob as glob
import pandas as pd
import numpy as np 

nlr_path=glob.glob('/global/scratch/users/chandlersutherland/e19/re_digest/*/nlr_seq_coverage.tsv')
gene_path=glob.glob('/global/scratch/users/chandlersutherland/e19/re_digest/*/gene_seq_coverage.tsv')

#CviAII_nlr=pd.read_csv(nlr_path[0], sep='\t', names=['chrom', 'start', 'stop', 'gene', 'cov'])
CviAII_nlr=pd.read_csv(nlr_path[0], sep='\t', names=['Chrom', 'start', 'stop', 'gene', 'strand', 'cov'])
CviAII_nlr['nlr_width']=CviAII_nlr['stop']-CviAII_nlr['start']
nlr_coverage=sum(CviAII_nlr['cov'])/sum(CviAII_nlr['nlr_width'])

def bed_coverage_calc(path):
    csv=pd.read_csv(path, sep='\t', names=['chrom', 'start', 'stop', 'gene', 'strand', 'cov'])
    csv['width']=csv['stop']-csv['start'] #calculate gene length 
    coverage=sum(csv['cov']/sum(csv['width'])) #calculate total % of bases covered 
    return coverage 

def coverage_calc(path):
    csv=pd.read_csv(path, sep='\t')

bed_result=[]
genome_result=pd.DataFrame()

for enzyme in ['AluI', 'BfaI', 'CviAII', 'CviQI', 'DpnII', 'FatI', 'HaeIII', 'HpyCH4V', 'MluCl', 'MseI', 'MspI', 'TaqI']:
    print(enzyme)
    nlr_path='/global/scratch/users/chandlersutherland/e19/re_digest/'+enzyme+'/nlr_seq_coverage.tsv'
    gene_path='/global/scratch/users/chandlersutherland/e19/re_digest/'+enzyme+'/gene_seq_coverage.tsv'
    genome_path='/global/scratch/users/chandlersutherland/e19/re_digest/'+enzyme+'/'+enzyme+'_seq_coverage.tsv'
    
    nlr_coverage=bed_coverage_calc(nlr_path)
    gene_coverage=bed_coverage_calc(gene_path)

    genome_coverage=pd.read_csv(genome_path, sep='\t')
    genome_coverage['enzyme']=enzyme 
    genome_result=pd.concat([genome_result,genome_coverage])
    result.append([enzyme, nlr_coverage, gene_coverage])

coverage_result=pd.DataFrame(result)
coverage_result.columns = ['enzyme', 'nlr_cov', 'gene_cov']
coverage_result.to_csv('/global/scratch/users/chandlersutherland/e19/re_digest/mapping_result.csv')