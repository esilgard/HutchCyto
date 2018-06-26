*__HutchCyto__*
==============================================================================

Primary directory for the Fred Hutch standard cytogenetics processing pipeline
------------------------------------------------------------------------------

* cyto_engine.py: main script - parse arguments and call other scripts
* command_line_flags.txt - keep track of required and optional command line arguments
	* -f: input file path
	* -o: output file path
	* -cl: clonality flag (t or f)
* global_strings.py: a set of string variables used in various scripts in this directory
* iscn_parser.py: a program to parse the cleaned karyotype string into a queryable data structure [ISCN formatting guidelines](http://www.cydas.org/Docs/ISCNAnalyser/Analysis.html)
* iscn_string_cleaner.py: a simple script to clean up the ISCN string by stripping off trailing or leading free text, accounting for common and unambiguous formatting irregularities
* parse_abnormalities.py: create a list of "types" of abnormalities by chromosome and arm (e.g. any deletion of 12p) with cell counts
* eln_classification.py: an ELN risk classification for AML based ONLY on the standard cytogenetics (real classifications should involve available FISH and molecular data as well)
* output_results.py: yeah...ya know...make it a pretty json (this is the best place to amend output format and/or location)
