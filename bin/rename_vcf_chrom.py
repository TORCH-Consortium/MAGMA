#! /usr/bin/env python

'''Original author Jody Phelan at https://github.com/jodyphelan/pathogen-profiler/blob/master/scripts/rename_vcf_chrom.py'''
import sys
import argparse
import subprocess

def errlog(x,ext=False):
    sys.stderr.write('\033[91m' + str(x) + '\033[0m' + '\n')
    if ext==True:
        quit(1)

def cmd_out(cmd,verbose=1):
    if verbose==2:
        sys.stderr.write("\nRunning command:\n%s\n" % cmd)
        stderr = open("/dev/stderr","w")
    elif verbose==1:
        sys.stderr.write("\nRunning command:\n%s\n" % cmd)
        stderr = open("/dev/null","w")
    else:
        stderr = open("/dev/null","w")
    try:
        res = subprocess.Popen(cmd,shell=True,stderr = stderr,stdout=subprocess.PIPE)
        for l in res.stdout:
            yield l.decode().rstrip()
    except:
        errlog("Command Failed! Please Check!")
        raise Exception
    stderr.close()

def main(args):
    generator = cmd_out("bcftools view " + args.vcf) if args.vcf else sys.stdin
    convert = dict(zip(args.source,args.target))

    # Open the output file for writing
    with open(args.outfile, 'w') as outfile:
        for l in generator:
            if l[0]=="#":
                outfile.write(l.strip()+"\n")
            else:
                row = l.strip().split()
                row[0] = convert[row[0]]
                outfile.write("\t".join(row)+"\n")

parser = argparse.ArgumentParser(description='tbprofiler script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--vcf',type=str,help='')
parser.add_argument('--source',nargs="+",type=str,help='')
parser.add_argument('--target',nargs="+",type=str,help='')
parser.add_argument('--outfile', type=str, required=True, help='Output VCF file')
parser.set_defaults(func=main)
args = parser.parse_args()
args.func(args)
