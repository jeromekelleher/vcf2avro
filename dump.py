"""
Short script to dump an Avro encoded VCF.
"""
from __future__ import print_function
from __future__ import division 

import sys

import avro
import avro.io
import avro.schema
import avro.datafile

def dump_file(f):
    reader = avro.datafile.DataFileReader(f, avro.io.DatumReader())
    for record in reader:
        print(record["POS"])

def main():
    if len(sys.argv) != 2:
        print("Usage: dump <avro_vcf>")
        sys.exit()
    with open(sys.argv[1]) as f:
        dump_file(f)

if __name__ == "__main__":
    main()
