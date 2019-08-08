''' author@esilgard '''
# Copyright (c) 2013-2018 Fred Hutchinson Cancer Research Center
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

import sys, os, csv
import output_results
from datetime import datetime
import global_strings as gb
import iscn_string_cleaner
import iscn_parser
import parse_abnormalities

'''
initial script of the standard cytogenetics engine do deal with command line parsing and module outputs
should exit with a non-zero status for any fatal errors and
output warnings and results in json format to CWD in the file provided in cmd line args
'''

## declare output dictionary for values, warnings, and metadata
OUTPUT_DICTIONARY = {}
OUTPUT_DICTIONARY[gb.ERRS] = []

## path to file containing command line flags and descriptions ##
## in the format -char<tab>description<tab>verbose_description(for help and error messages) ##
try:
    COMMAND_LINE_FLAG_FILE = open('command_line_flags.txt', 'r')
    try:
        ## set of required flags for program to run successfully ##
        REQUIRED_FLAGS = set([])
        ## dictionary of actual flags:argument values ##
        ARGUMENTS = {}
        ## dictionary of flag:tuple(flag description,verbose flag description) ##
        COMMAND_LINE_FLAGS = {}
        for line in COMMAND_LINE_FLAG_FILE.readlines():
            line = line.strip().split('\t')
            if line[1] == 'required':
                REQUIRED_FLAGS.add(line[0])
            COMMAND_LINE_FLAGS[line[0]] = (line[2], line[3])
        COMMAND_LINE_FLAG_FILE.close()
        ARGS = sys.argv[1:]
    except IOError:
        sys.stderr.write('FATAL ERROR: command line flag dictionary could not be established \
                        from file, potential formatting error.  program aborted.')
        sys.exit(1)
except EnvironmentError:
    sys.stderr.write('FATAL ERROR: command line flag file not found.  program aborted.')
    sys.exit(1)

## parse the ARGUMENTS from arg1 on into a dictionary - notify user of unrecognized flags
## NOTE - this does assume that flags start in the first position
## and every other argument is a flag
for index in range(0, len(ARGS)-1, 2):
    if ARGS[index] in COMMAND_LINE_FLAGS:
        ARGUMENTS[ARGS[index]] = ARGS[index+1]
    else:
        OUTPUT_DICTIONARY[gb.ERRS].append({gb.ERR_TYPE: 'Warning', gb.ERR_STR: 'nonfatal error: \
        unrecognized flag: ' + ARGS[index] + ', this flag will not be excluded. Refer to \
        command_line_flags.txt for a complete list and description of command line flags'})

## build the dictionary for the json output ##
OUTPUT_DICTIONARY[gb.CNTL] = {}
OUTPUT_DICTIONARY[gb.CNTL][gb.DATE] = str(datetime.today().isoformat())
OUTPUT_DICTIONARY[gb.REPORTS] = []
OUTPUT_DICTIONARY[gb.CNTL][gb.INPUT_F] = ARGUMENTS.get('-f')
OUTPUT_DICTIONARY[gb.CNTL][gb.CLONE_FLG] = True
if ARGUMENTS.get('-cl') == 'n':
    OUTPUT_DICTIONARY[gb.CNTL][gb.CLONE_FLG] = False

## ERR out for missing flags that are required ##
MISSING_FLAGS = REQUIRED_FLAGS-set(ARGUMENTS.keys())
if len(MISSING_FLAGS) > 0:
    for each_flag in MISSING_FLAGS:
        sys.stderr.write('FATAL ERROR: missing required flag: ' + each_flag + ' ' + COMMAND_LINE_FLAGS[each_flag][1] + os.linesep)
    sys.exit(1)
else:
    with open(ARGUMENTS.get('-f'), mode='r') as infile:
        reader = csv.DictReader(infile) #, delimiter='\t')
        raw_cyto_d = dict((rows[gb.KARYO_ID],rows[gb.KARYO_STR]) for rows in reader)
    for k_id,k_str in raw_cyto_d.items():
        clean_str = iscn_string_cleaner.get(k_str)
        cell_list = iscn_parser.get(clean_str)

        OUTPUT_DICTIONARY[gb.REPORTS].append({k_id:parse_abnormalities.get(cell_list, clean_str, OUTPUT_DICTIONARY[gb.CNTL][gb.CLONE_FLG])})
    ## output results to file ##
    OUTPUT_RETURN = output_results.main(ARGUMENTS.get('-o'), OUTPUT_DICTIONARY)
    if OUTPUT_RETURN:
        sys.exit(1)

