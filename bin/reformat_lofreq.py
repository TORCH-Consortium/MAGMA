import ast
import argparse

import pandas as pd

#NOTE: Comment out the multi-line string as it's raising a runtime error and isn't used in the script

# header_formats = """##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
# ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
# ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
# ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
# ##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Minimum DP observed within the GVCF block">
# ##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
# """

def read_vcf(filename):
    with open(filename, 'r') as vcf:
        header = ''
        for line in vcf:
            if line[:2] != '##':
                break
            else:
                header += line
        df = pd.read_csv(vcf, header=None, sep='\t')
    df.columns = line[:-1].split('\t')
    return df, header

def write_vcf(filename, df, header):
    with open(filename, 'w') as vcf:
        vcf.write(header)
        df.to_csv(vcf, sep='\t', index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Analyse resistance output from XBS Pipeline')
    parser.add_argument('lofreq_vcf_file', metavar='lofreq_vcf_file', type=str, help='The input lofreq vcf file')
    parser.add_argument('lofreq_sample_name', metavar='lofreq_sample_name', type=str, help='The sample name')
    parser.add_argument('outfile', metavar='outfile', type=str, help='The name of the output VCF file')

    args = vars(parser.parse_args())

    vcf, header = read_vcf(args['lofreq_vcf_file'])
    vcf['FORMAT'] = 'GT:AD:DP:GQ:PL'

    for idx, row in vcf.iterrows():
        info = [ast.literal_eval(i.split('=')[1]) for i in row['INFO'].split(';')[:4]]
        ref_dp = sum(info[3][:2])
        alt_dp = sum(info[3][2:])
        GT = 1
        AD = '{},{}'.format(ref_dp, alt_dp)
        DP = sum(info[3])
        GQ = 99
        PL = '1800,0'
        vcf.loc[idx, args['lofreq_sample_name']] = '{}:{}:{}:{}:{}'.format(GT,AD,DP,GQ,PL)

    write_vcf(args['outfile'], vcf, header)
