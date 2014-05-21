/*
** Copyright (C) 2014, Jerome Kelleher 
**  
** This file is part of vcf2avro.
** 
** Wormtable is free software: you can redistribute it and/or modify
** it under the terms of the GNU Lesser General Public License as published by
** the Free Software Foundation, either version 3 of the License, or
** (at your option) any later version.
** 
** Wormtable is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU Lesser General Public License for more details.
** 
** You should have received a copy of the GNU Lesser General Public License
** along with vcf2avro.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <avro.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

static void 
fatal_error(char *msg, ...)
{
    va_list argp;
    fprintf(stderr, "vcfcat:");
    va_start(argp, msg);
    vfprintf(stderr, msg, argp);
    va_end(argp);
    fprintf(stderr, "\n");
    exit(EXIT_FAILURE);
}



/* 
 * Reads the specified column from the specified VCF file in 
 * avro format. This method uses the legacy datum API.
 */
void
read_data_legacy(const char *vcffile, const char *column)
{
    avro_file_reader_t file_reader;
    	
    printf("reading %s from %s\n", column, vcffile); 
    if (avro_file_reader(vcffile, &file_reader)) {
        fatal_error("Error opening file '%s': %s", vcffile, avro_strerror());
    }

	avro_file_reader_close(file_reader);
}

int 
main(int argc, char **argv)
{
	const char *vcffile, *column;
    if (argc != 3) {
		fprintf(stderr, "usage: %s <FILENAME> <COLUMN>\n", argv[0]);
		exit(EXIT_FAILURE);
    }
    vcffile = argv[1];
    column = argv[2];
    read_data_legacy(vcffile, column);
    return EXIT_SUCCESS;
}
