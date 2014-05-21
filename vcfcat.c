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

/* for convenience until we change the avro type to string */
#define MAX_STRING 8192

/* 
 * Prints out the values in the specified row.
 */
static void
print_row_legacy(avro_datum_t row, const char *col)
{
    int32_t i32;
    int64_t i64;
    float  f;
    char *str;
    char dest[MAX_STRING];
    avro_datum_t union_dt, value_dt;

    if (avro_record_get(row, col, &union_dt) != 0) {
        fatal_error("Cannot get union from row: %s", avro_strerror());
    }
    if (is_avro_union(union_dt)) {
        value_dt = avro_union_current_branch(union_dt);
        if (value_dt == NULL) {
            fatal_error("Error getting union branch");
        }
    } else {
        fatal_error("case not handled");
    }
    if (is_avro_int32(value_dt)) {
        if (avro_int32_get(value_dt, &i32) != 0) {
            fatal_error("Error reading int: %s", avro_strerror());
        }
        printf("%d", i32);
    } else if (is_avro_bytes(value_dt)) {
        if (avro_bytes_get(value_dt, &str, &i64) != 0) {
            fatal_error("Error converting to bytes");
        }
        if (i64 >= MAX_STRING) {
            fatal_error("max string size exceeded");
        }
        memcpy(dest, str, i64);
        dest[i64] = '\0';
        printf("%s", dest);
    } else if (is_avro_float(value_dt)) {
        if (avro_float_get(value_dt, &f) != 0) {
            fatal_error("Error converting to float");
        }
        printf("%f", f);
    } else if (is_avro_null(value_dt)) {
        printf("NA");
    } else {
        fatal_error("Avro type not handled");
    }
    printf("\n"); 
}

/* 
 * Reads the specified column from the specified VCF file in 
 * avro format. This method uses the legacy datum API.
 */
static void
read_data_legacy(const char *vcffile, const char *column)
{
    int ret = 0;
    avro_file_reader_t file_reader;
    avro_schema_t read_schema, subset_schema, col_schema;
    avro_datum_t row;

    if (avro_file_reader(vcffile, &file_reader)) {
        fatal_error("Error opening file '%s': %s", vcffile, avro_strerror());
    }
    read_schema = avro_file_reader_get_writer_schema(file_reader); 
    if (read_schema == NULL) {
        fatal_error("Error getting read_schema: %s", avro_strerror());
    }
    col_schema = avro_schema_get_subschema(read_schema, column); 
    if (col_schema == NULL) {
        fatal_error("Error getting column '%s': %s", column, avro_strerror());
    }
    /* Create the subset schema */
    subset_schema = avro_schema_record("VCF", "vcf.avro");
    if (subset_schema == NULL) {
        fatal_error("Error making subset schema:%s", avro_strerror());
    }
    if (avro_schema_record_field_append(subset_schema, column, col_schema) 
            != 0) {
        fatal_error("Error building subset schema:%s", avro_strerror());
    }
	ret = avro_file_reader_read(file_reader, subset_schema, &row);
    while (ret == 0) {
        print_row_legacy(row, column);
        avro_datum_decref(row);
        ret = avro_file_reader_read(file_reader, subset_schema, &row);
    }
    /* For some reason, the last call to reader_read returns 28, when 
     * it should be returning EOF. No errors appeared to have happened,
     * so we just work around this here for now
     */
    if (ret != EOF && ret != 28) {
        fatal_error("Error reading records: %s", avro_strerror());
    }
    avro_schema_decref(subset_schema);
    avro_schema_decref(read_schema);
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
