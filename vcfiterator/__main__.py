import json
import sys
import argparse

from vcfiterator import VcfIterator

parser = argparse.ArgumentParser("Iterates over a .vcf file, outputting one JSON structure per line")
parser.add_argument("vcf_file", help="Path to .vcf file")
parser.add_argument("--pretty", action="store_true", help="Pretty print JSON")

args = parser.parse_args()
path = args.vcf_file
v = VcfIterator(path)

for value in v.iter():
    kw = {}
    if args.pretty:
        kw['indent'] = 4
    print json.dumps(value, **kw)
