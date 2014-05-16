#!python
# 
# Copyright (C) 2014, Jerome Kelleher jk@well.ox.ac.uk 
#
# This file is part of vcf2avro.
# 
# Wormtable is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# Wormtable is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public License
# along with vcf2avro.  If not, see <http://www.gnu.org/licenses/>.
# 
"""
VCF processing for wormtable. 

Implementation Note: We use bytes throughout the parsing process here for
a few reasons. Mostly, this is because it's much easier to deal with bytes
values within the C module, as we'd have to decode Unicode objects before 
getting string. At the same time, it's probably quite a bit more efficient 
to work with bytes directly, so we win both ways. It's a bit tedious making 
sure that all the literals have a 'b' in front of them, but worth the 
effort.
"""
from __future__ import print_function
from __future__ import division 

import os
import sys
import time
import shutil 
import argparse
import tempfile

import avro
import avro.io
import avro.schema
import avro.datafile

# VCF Fixed columns

CHROM_NAME = b"CHROM"
POS_NAME = b"POS"
ID_NAME = b"ID"
REF_NAME = b"REF"
ALT_NAME = b"ALT"
QUAL_NAME = b"QUAL"
FILTER_NAME = b"FILTER"
INFO_NAME = b"INFO"

VCF_FIXED_COLUMNS = [CHROM_NAME, POS_NAME, ID_NAME, REF_NAME, ALT_NAME, 
        QUAL_NAME, FILTER_NAME]

# Descriptions taken from 
# http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41 
CHROM_DESCRIPTION = b"""chromosome: an identifier from the reference genome or an \
angle-bracketed ID String ("<ID>") pointing to a contig in the assembly file"""
POS_DESCRIPTION = b"""position: The reference position, with the 1st base having \
position 1"""
ID_DESCRIPTION = b"""semi-colon separated list of unique identifiers where \
available"""
REF_DESCRIPTION = b"""reference base(s): Each base must be one of A,C,G,T,N \
(case insensitive)""" 
ALT_DESCRIPTION = b"""comma separated list of alternate non-reference alleles\
called on at least one of the samples"""
QUAL_DESCRIPTION = b"""phred-scaled quality score for the assertion made in \
ALT. i.e. -10log_10 prob(call in ALT is wrong)."""
FILTER_DESCRIPTION = b"""PASS if this position has passed all filters, i.e. \
a call is made at this position. Otherwise, if the site has not passed all \
filters, a semicolon-separated list of codes for filters that fail. """

# char used to seperate VCF columns from their prefix, e.g INFO.AF
COLUMN_SEPARATOR = b"_"

# Special values in VCF
MISSING_VALUE = b"."

# Strings used in the header for identifiers
ID = b"ID"
INFO = b"INFO"
DESCRIPTION = b"Description"
NUMBER = b"Number"
TYPE = b"Type"
INTEGER = b"Integer"
FLOAT = b"Float"
FLAG = b"Flag"
CHARACTER = b"Character"
STRING = b"String"

VARIABLE_SIZE = 0 

class ProgressMonitor(object):
    """
    Class representing a progress monitor for a terminal based interface.
    """
    def __init__(self, total, units):
        self.__total = total
        self.__units = units
        self.__progress_width = 40
        self.__bar_index = 0
        self.__bars = "/-\\|"
        self.__start_time = time.clock()

    def update(self, processed):
        """
        Updates this progress monitor to display the specified number 
        of processed items.
        """
        complete = processed / self.__total
        filled = int(complete * self.__progress_width)
        spaces = self.__progress_width - filled 
        bar = self.__bars[self.__bar_index]
        self.__bar_index = (self.__bar_index + 1) % len(self.__bars)
        elapsed = max(1, time.clock() - self.__start_time)
        rate = processed / elapsed
        s = '\r[{0}{1}] {2:5.1f}% @{3:8.1E} {4}/s {5}'.format('#' * filled, 
            ' ' * spaces, complete * 100, rate, self.__units, bar)
        sys.stdout.write(s)
        sys.stdout.flush()
         
    def finish(self):
        """
        Completes the progress monitor.
        """
        print()



class FileReader(object):
    """
    A class for reading data files from a variety of sources and 
    with progress updating.
    """
    def __init__(self, in_file):
        if in_file == '-':
            self.__input_file = sys.stdin
            if sys.version_info[:2] >= (3, 1):
                try:
                    self.__input_file = sys.stdin.buffer
                except AttributeError:
                    # When we're testing, we replace stdin with an ordinary
                    # file. This does not support buffer, but is also not 
                    # needed so we can skip this step
                    pass
            self.__progress_monitor = None
            self.__input_file_size = None
            self.__progress_file = None
        else:
            if in_file.endswith(".gz"):
                # Detect broken GZIP handling in 2.7/3.2 and others and abort
                # TODO this has been fixed upstream and can be removed at 
                # some point.
                f = gzip.open(in_file, "rb")
                try:
                    s = f.readline()
                except Exception as e:
                    print(BROKEN_GZIP_MESSAGE)
                    sys.exit(1)
                f.close() 
                # Carry on as before
                self.__input_file = gzip.open(in_file, "rb")
                self.__progress_file = self.__input_file.fileobj
            else:
                self.__input_file = open(in_file, "rb")
                self.__progress_file = self.__input_file 
            statinfo = os.stat(in_file)
            self.__input_file_size = statinfo.st_size 
        self.__progress_update_rows = 2**32 
        self.__progress_monitor = None

    def get_progress_update_rows(self):
        """
        Returns the number of rows after which we should update the progress 
        monitor.
        """
        return self.__progress_update_rows

    def set_progress_update_rows(self, update_rows):
        """
        Sets the number of rows after which we should update progress 
        to the specified value.
        """
        self.__progress_update_rows = update_rows

    def get_input_file(self):
        """
        Returns the File object that is the source of the information 
        in this reader.
        """
        return self.__input_file

    def set_progress(self, progress):
        """
        If progress is True turn on progress monitoring for this GTF reader.
        """
        if progress:
            self.__progress_monitor = ProgressMonitor(self.__input_file_size, 
                    "bytes")
            self.__progress_monitor.update(0)
            self.__progress_update_rows = 100
            if self.__input_file_size > 2**30:
                self.__progress_update_rows = 1000
    
    def update_progress(self):
        """
        Reads the position we are at in the underlying file and uses this to 
        update the progress bar, if used.
        """
        if self.__progress_monitor is not None:
            t = self.__progress_file.tell() 
            self.__progress_monitor.update(t)

    def finish_progress(self):
        """
        Finishes up the progress monitor, if in use.
        """ 
        if self.__progress_monitor is not None:
            self.update_progress()
            self.__progress_monitor.finish()

    def close(self):
        """
        Closes any open files on this Reader.
        """
        self.__input_file.close()
        if self.__progress_file is not None: 
            self.__progress_file.close()

 

class VCFReader(FileReader):
    """
    A class for reading VCF files. 
    """
    def __init__(self, vcf_file):
        super(VCFReader, self).__init__(vcf_file)
        self.__genotypes = []
        self.__truncate = False
        self.read_header()

    def set_truncate_REF_ALT(self, truncate):
        """
        If true, truncate REF and ALT columns to be no more than 253 characters 
        long. This is a temporary workaround until more sophisticated truncation 
        across all columns is implemented.
        """
        self.__truncate = truncate

    def parse_version(self, s):
        """
        Parse the VCF version number from the specified string.
        """
        self._version = -1.0
        tokens = s.split(b"v")
        if len(tokens) == 2:
            self._version = float(tokens[1])

    def parse_header_line(self, s):
        """
        Processes the specified header string to get the genotype labels.
        """
        self.__genotypes = s.split()[9:]

    def _get_converter(self, avro_type, num_elements):
        """
        Returns a conversion function for the specified type and number 
        of elements.
        """
        d = {"boolean":int, "int": int, "float": float, "bytes": None}
        if num_elements == 1:
            f = d[avro_type]
        else:
            g = d[avro_type] 
            def conv(s):
                ret = []
                for tok in s.split(","):
                    ret.append(g(tok))
                return ret
            if g is None:
                f = None
            else:
                f = conv
        return f

    def add_column_definition(self, name, description, avro_type, num_elements=1):
        s = """{{"name": "{0}", """.format(name)
        if num_elements == 1 or avro_type == "bytes":
            t = "\"{0}\"".format(avro_type)
        else:
            t = """{{"type":"array", "items":"{0}"}}""".format(avro_type)
        s += """"type": [{0}, "null"]}}, """.format(t)
        f = self._get_converter(avro_type, num_elements)
        self.__columns[name] = f 
        self.__schema += s + "\n"

    def add_int_column(self, name, description):
        """
        Adds an int column to the Avro schema.
        """
        self.add_column_definition(name, description, "int")

    def add_char_column(self, name, description):
        """
        Adds a char column to the Avro schema.
        """
        self.add_column_definition(name, description, "bytes")
    
    def add_float_column(self, name, description):
        """
        Adds a float column to the Avro schema.
        """
        self.add_column_definition(name, description, "float")


    def add_column(self, prefix, line):
        """
        Adds a VCF column using the specified metadata line with the specified 
        name prefix to the specified table.
        """
        d = {}
        s = line[line.find(b"<") + 1: line.find(b">")]
        for j in range(3):
            k = s.find(b",")
            tokens = s[:k].split(b"=")
            s = s[k + 1:]
            d[tokens[0]] = tokens[1]
        tokens = s.split(b"=", 1)
        d[tokens[0]] = tokens[1]
        name = d[ID]
        description = d[DESCRIPTION].strip(b"\"")
        number = d[NUMBER]
        num_elements = VARIABLE_SIZE 
        try:
            # If we can parse it into a number, do so. If this fails than use
            # a variable number of elements.
            num_elements = int(number)    
        except ValueError as v:
            pass
        # We can also have negative num_elements to indicate variable column
        if num_elements < 0:
            num_elements = VARIABLE_SIZE 
        st = d[TYPE]
        if st == INTEGER:
            element_type = "int" 
            element_size = 4
        elif st == FLOAT: 
            element_type = "float"
            element_size = 4
        elif st == FLAG: 
            element_type = "int" # TODO change to boolean
            element_size = 1
            num_elements = 1
        elif st == CHARACTER: 
            element_type = "bytes" 
            element_size = 1
        elif st == STRING: 
            num_elements = VARIABLE_SIZE 
            element_type = "bytes" 
            element_size = 1
        else:
            raise ValueError("Unknown VCF type:", st)
        
        self.add_column_definition(prefix + COLUMN_SEPARATOR + name, 
                description, element_type, num_elements)

    def generate_schema(self):
        """
        Reads the header from the specified VCF file and returns an 
        avro schema in JSON format.
        """
        self.__schema = """
        {"namespace":"vcf.avro",
         "type":"record",
         "name":"VCF",
         "fields": [
        """
        self.__columns = {}
        info_descriptions = []
        genotype_descriptions = []
        
        if self._version < 4.0:
            raise ValueError("VCF versions < 4.0 not supported")
        for s in self.__header:
            # skip FILTER values 
            if s.startswith(b"##INFO"):
                info_descriptions.append(s)
            elif s.startswith(b"##FORMAT"):
                genotype_descriptions.append(s)

        # Add the fixed columns
        self.add_char_column(CHROM_NAME, CHROM_DESCRIPTION)
        self.add_int_column(POS_NAME, POS_DESCRIPTION)
        self.add_char_column(ID_NAME, ID_DESCRIPTION)
        self.add_char_column(REF_NAME, REF_DESCRIPTION)
        self.add_char_column(ALT_NAME, ALT_DESCRIPTION)
        self.add_float_column(QUAL_NAME, QUAL_DESCRIPTION)
        self.add_char_column(FILTER_NAME, FILTER_DESCRIPTION)
        
        for s in info_descriptions:
            self.add_column(INFO_NAME, s)
        for genotype in self.__genotypes:
            for s in genotype_descriptions: 
                self.add_column(genotype, s) 
        # Remove the last ','
        self.__schema = self.__schema.rstrip("\n ,")
        self.__schema += "]}"
        return self.__schema, self.__columns

    def read_header(self):
        """
        Read header lines, parse version and column names
        """
        f = self.get_input_file()
        self.__header = [f.readline()]  
        while self.__header[-1].startswith(b"##"):
            self.__header.append(f.readline())
        self.parse_version(self.__header[0])
        self.parse_header_line(self.__header.pop()) 
        

    def rows(self, table_columns):
        """
        Returns an iterator over the rows in this VCF file. Each row is a 
        dictionary mapping column positions to their encoded string values.
        """
        # First we construct the mappings from the various parts of the 
        # VCF row to the corresponding column name in the Avro schema 
        num_columns = len(table_columns) 
        all_fixed_columns = VCF_FIXED_COLUMNS 
        fixed_columns = []
        # weed out the columns that are not in the table
        for j in range(len(all_fixed_columns)):
            name = all_fixed_columns[j]
            if name in table_columns:
                fixed_columns.append((j, name))
        info_columns = {}
        genotype_columns = [{} for g in self.__genotypes]
        for k in table_columns: 
            if COLUMN_SEPARATOR in k:
                split = k.split(COLUMN_SEPARATOR)
                if split[0] == INFO:
                    name = COLUMN_SEPARATOR.join(split[1:])
                    info_columns[name] = k 
                else:
                    g = COLUMN_SEPARATOR.join(split[:-1])
                    name = split[-1]
                    index = self.__genotypes.index(g)
                    genotype_columns[index][name] = k 
        ref_index = 3
        alt_index = 4
        # Now we are ready to process the file.
        update_rows = self.get_progress_update_rows()
        num_rows = 0
        for s in self.get_input_file():
            row = {} 
            l = s.split()
            # Read in the fixed columns
            for vcf_index, wt_index in fixed_columns:
                if l[vcf_index] != MISSING_VALUE:
                    row[wt_index] = l[vcf_index]
            # Now process the info columns.
            for mapping in l[7].split(b";"):
                tokens = mapping.split(b"=")
                name = tokens[0]
                if name in info_columns:
                    col = info_columns[name]
                    if len(tokens) == 2:
                        row[col] = tokens[1]
                    else:
                        # This is a Flag column.
                        row[col] = b"1"
            # Process the genotype columns, if they exist
            if len(l) > 8:
                j = 0
                fmt = l[8].split(b":")
                for genotype_values in l[9:]:
                    tokens = genotype_values.split(b":")
                    if len(tokens) == len(fmt):
                        for k in range(len(fmt)):
                            if fmt[k] in genotype_columns[j]:
                                col = genotype_columns[j][fmt[k]]
                                tok = tokens[k]
                                # FIXME this is a hack to detect missing values 
                                # in genotype columns. I'm not sure why anybody 
                                # would do this, but we need it to parse the 
                                # example VCF from the 1000genomes site.
                                if tok != MISSING_VALUE and tok != b".,.":
                                    row[col] = tok 
                    j += 1
            d = {}
            """
            for k  in table_columns:
                if k not in row:
                    print("Column", k, "not found")
                    row[k] = "0"
            """
            for k, v in row.items():
                f = table_columns[k]
                if f is None:
                    d[k] = v
                else:
                    d[k] = f(v)
            yield d
            num_rows += 1
            if num_rows % update_rows == 0:
                self.update_progress()
        self.finish_progress()


class ProgramRunner(object):
    """
    Class responsible for running the vcf2wt program.
    """
    def __init__(self, args):
        self.__destination = args.DEST 
        self.__db_cache_size = args.cache_size
        self.__force = args.force
        self.__generate_schema = args.generate_schema
        self.__progress = not args.quiet
        self.__quiet = args.quiet
        self.__schema = args.schema
        self.__truncate = args.truncate
        self.__tmp_dirs = []
        self.__tmp_files = []
        self.__table = None
        self.__column_map = None
        self.__reader = VCFReader(args.SOURCE)
        self.__writer = None
        # if reading from STDIN, set progress monitor to False regardless
        if args.SOURCE == '-': 
            self.__progress = False

       
    def get_table(self):
        return self.__table

    def generate_schema(self):
        """
        Reads the header of the input VCF and generates a schema file.
        """
        self.__schema, columns = self.__reader.generate_schema()
        self.__column_map = columns
        #with open("schema.json", "w") as f:
        #    f.write(self.__schema)
        #with open("schema.json", "r") as f:
        #    self.__schema = f.read()

    def create_table(self):
        """
        Creates the table and reads the column information for the VCF reader.
        """
        schema = avro.schema.parse(self.__schema)
        self.__output_file = open(self.__destination, "w")
        self.__writer = avro.datafile.DataFileWriter(self.__output_file, 
                avro.io.DatumWriter(), schema, codec="deflate")

    def write_table(self):
        """
        Writes the table, assuming that we have created a directory with
        a table ready for writing.
        """
        self.__reader.set_progress(self.__progress)
        self.__reader.set_truncate_REF_ALT(self.__truncate)
        for r in self.__reader.rows(self.__column_map):
            self.__writer.append(r) 
        self.__reader.close()
        self.__reader = None
        self.__writer.close()
        self.__writer = None

    def run(self):
        """
        Top level entry point.
        """ 
        if self.__schema is None:
            self.generate_schema()
        
        if os.path.exists(self.__destination):
            if self.__force:
                os.unlink(self.__destination)
            else:
                s = "'{0}' exists; use -f to overwrite".format(self.__destination)
                self.error(s)
                
        if self.__generate_schema:
            # copy the schema and we're done.
            shutil.copyfile(self.__schema, self.__destination) 
        else:
            self.create_table()
            self.write_table()

    def error(self, s):
        """
        Raises and error and exits.
        """
        print("Error:", s)
        sys.exit(0)

    def cleanup(self):
        """
        Cleans up any temporary files shuts down any running processes.
        """
        for f in self.__tmp_dirs:
            shutil.rmtree(f)
        for f in self.__tmp_files:
            os.unlink(f)
        if self.__reader is not None:
            self.__reader.close()
        if self.__writer is not None:
            self.__writer.close()

def main(args=None):
    prog_description = "Convert a VCF file to Wormtable format."
    parser = argparse.ArgumentParser(description=prog_description) 
    parser.add_argument("SOURCE", 
        help="VCF file to convert (use '-' for STDIN)")   
    parser.add_argument("DEST", 
        help="""Output wormtable home directory, or schema file 
            if we are generating a candidate schema using the 
            --generate-schema option""")   
    parser.add_argument("--quiet", "-q", action="store_true", 
        default=False,
        help="Suppress progress monitor")   
    parser.add_argument("--force", "-f", action="store_true", default=False,
        help="Force over-writing of existing wormtable")   
    parser.add_argument("--truncate", "-t", action="store_true", default=False,
        help="""Truncate values that are too large for a column and store
            the maximum the column will allow. Currently this option 
            only supports truncating REF and ALT column values more than 254 
            characters long. REF and ALT values are truncated to 253 characters
            and suffixed with a '+' to indicate that truncation has 
            occured""")   
    parser.add_argument("--cache-size", "-c", default="64M",
        help="cache size in bytes; suffixes K, M and G also supported.")   
    g = parser.add_mutually_exclusive_group()
    g.add_argument("--generate-schema", "-g", action="store_true", 
        default=False,
        help="""Generate a schema for the source VCF file and 
            write to DEST. Only reads the header of the VCF file.""")   
    g.add_argument("--schema", "-s", default=None,
        help="""Use schema from the file SCHEMA rather than default 
                generated schema""")
    parsed_args = parser.parse_args(args)
    runner = ProgramRunner(parsed_args) 
    try:
        runner.run()
    finally:
        runner.cleanup()

if __name__ == "__main__":
    main()

