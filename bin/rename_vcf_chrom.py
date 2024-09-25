#! /usr/bin/env python
#
# Copyright (c) 2021-2024 MAGMA pipeline authors, see https://doi.org/10.1371/journal.pcbi.1011648
#
# This file is part of MAGMA pipeline, see https://github.com/TORCH-Consortium/MAGMA
#
# For quick overview of GPL-3 license, please refer
# https://www.tldrlegal.com/license/gnu-general-public-license-v3-gpl-3
#
# - You MUST keep this license with original authors in your copy
# - You MUST acknowledge the original source of this software
# - You MUST state significant changes made to the original software
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program . If not, see <http://www.gnu.org/licenses/>.
#

'''Original author Jody Phelan at https://github.com/jodyphelan/pathogen-profiler/blob/master/scripts/rename_vcf_chrom.py'''
import sys
import argparse
import subprocess

def errlog(x,ext=False):
    sys.stderr.write('\033[91m' + str(x) + '\033[0m' + '\n')
    if ext==True:
        quit(1)

def cmd_out(cmd,verbose=1):
    cmd = "set -u pipefail; " + cmd
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
    for l in generator:
        if l[0]=="#":
            sys.stdout.write(l.strip()+"\n")
        else:
            row = l.strip().split()
            row[0] = convert[row[0]]
            sys.stdout.write("\t".join(row)+"\n")

parser = argparse.ArgumentParser(description='tbprofiler script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--vcf',type=str,help='')
parser.add_argument('--source',nargs="+",type=str,help='')
parser.add_argument('--target',nargs="+",type=str,help='')
parser.set_defaults(func=main)
args = parser.parse_args()
args.func(args)
