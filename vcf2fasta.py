#!/home/rs2474/anaconda3/envs/imktData/bin/python

### Convert VCF ouput from SLiM to fasta file for iMKT analysis

import pandas as pd
import numpy as np
import sys


def output_ref_sub_seq(slim_out, gen):
    with open(slim_out, 'r') as f:
        f_lines = f.readlines()
    ## get the original reference of bam CDS
    ref_line = list(filter(lambda line: line.startswith('Ancestral:'), f_lines))[0]
    ref_bam_cds = ref_line.split(':')[1].strip()
    ref_bam_cds_codons = [ref_bam_cds[i:i+3] for i in range(0, len(ref_bam_cds), 3)]
    
    ## check whether new reference seq exists
    new_ref_line = list(filter(lambda line: line.startswith('Generation {} New Ancestral: '.format(gen)), f_lines))[0]
    new_ref_cds = new_ref_line.split(':')[1].strip()
    new_ref_bam_cds_codons = [new_ref_cds[i:i+3] for i in range(0, len(new_ref_cds), 3)]
    ## check the number of complex codons in the simulation
    codon_diff = list(map(lambda ref, sub: sum(np.char.equal(np.array(list(ref)), np.array(list(sub)))), ref_bam_cds_codons,new_ref_bam_cds_codons))
    print("{} complex codons detected in counting divergence.".format(len(list(filter(lambda x: x < 2, codon_diff)))))
    return ref_bam_cds, new_ref_cds
    
    
def output_haplo_from_vcf(ref_seq, vcf_file):
    new_ref_cds = list(ref_seq)

    vcf_df = pd.read_csv(vcf_file, sep='\t', skiprows=14, header=0)
    seq_dict_hap1 = {k+'-1':new_ref_cds.copy() for k in vcf_df.iloc[:,9:].columns}
    seq_dict_hap2 = {k+'-2':new_ref_cds.copy() for k in vcf_df.iloc[:,9:].columns}
    
    ref_codon_pos = [list(range(i, i+3)) for i in range(0, len(ref_seq), 3)]
    snp_pos = (vcf_df['POS']-1).tolist() # VCF 1-index while SLiM uses 0-index
    
    ## check the number of complex codons in the simulation
    snp_per_codon = list(map(lambda codon_pos: len(set(snp_pos).intersection(codon_pos)), ref_codon_pos))
    print("{} potential complex codons detected in counting polymorphisms.".format(len(list(filter(lambda x: x >= 2, snp_per_codon)))))
    
    for row_idx, row in vcf_df.iterrows():
        idx = row['POS']-1
        ## check if the ref match up
        ## if the assertion fails, check whether you replace the substitutions in the ancestral sequence
        assert row['REF'] == new_ref_cds[idx]
        ## get the alt alleles
        alt_alleles = [row['ALT']] if (len(row['ALT'])==1) else row['ALT'].split(',')
        alt_allele_dict = dict(zip(range(1,len(alt_alleles)+1), alt_alleles))
        alt_allele_dict[0] = row['REF']

        for hap in vcf_df.iloc[:,9:]:
            seq_dict_hap1[hap+'-1'][idx] = alt_allele_dict[int(row[hap].split('|')[0])]
            seq_dict_hap2[hap+'-2'][idx] = alt_allele_dict[int(row[hap].split('|')[1])]

    seq_dict_hap1 = list(map(lambda x: (x[0],''.join(x[1])), seq_dict_hap1.items()))
    seq_dict_hap2 = list(map(lambda x: (x[0],''.join(x[1])), seq_dict_hap2.items()))
    
    return seq_dict_hap1, seq_dict_hap2


def output2fasta(fasta_file, new_ref_cds, ref_bam_cds, seq_dict_hap1, seq_dict_hap2):
    with open(fasta_file, 'w') as f:
        f.write('>REF\n'+new_ref_cds+'\n')
        for hap1, hap2 in zip(seq_dict_hap1, seq_dict_hap2):
            f.write('>'+hap1[0]+'\n'+hap1[1]+'\n')
            f.write('>'+hap2[0]+'\n'+hap2[1]+'\n')
        f.write('>OUTGROUP\n'+ref_bam_cds+'\n')
        

def main():
    vcf_file = sys.argv[1]
    slim_out = sys.argv[2]
    gen = int(list(filter(lambda x: 'Gen' in x, '.'.join(vcf_file.split('.')[:-1]).split('_')))[0][3:])
    fasta_file = '.'.join(vcf_file.split('.')[:-1]) + '.fa'
    
    print()
    print('From VCF file:', vcf_file)
    print('From SLiM output:', slim_out)
    
    ref_bam_cds, new_ref_cds = output_ref_sub_seq(slim_out, int(gen))
    seq_dict_hap1, seq_dict_hap2 = output_haplo_from_vcf(new_ref_cds, vcf_file)
    output2fasta(fasta_file, new_ref_cds, ref_bam_cds, seq_dict_hap1, seq_dict_hap2)

if __name__ == "__main__":
    main()
