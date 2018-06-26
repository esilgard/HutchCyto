'''author@esilgard'''
#
# Copyright (c) 2015-2018 Fred Hutchinson Cancer Research Center
#
# Licensed under the Apache License, Version 2.0: http://www.apache.org/licenses/LICENSE-2.0
#
import re
import global_strings as gb

__version__ = 'iscn_string_cleaner1.0'

def get(karyo_string):
    '''
    clean up the ISCN string that contains the ISCN karyotype info
    look for common formatting inconsistenies or typos in the karyotype
    return the cleaned string
    '''

    ## cut off FISH results that may have different formatting and vocabulary
    if 'nuc' in karyo_string:
        karyo_string = karyo_string[:karyo_string.find('nuc')]
    if 'NUC' in karyo_string:
        karyo_string = karyo_string[:karyo_string.find('NUC')]
        
    ## account for strings with ":" (instead of ";") only in the context of digit;digit
    if re.search(r'[\d]:[\d]', karyo_string):
        typo = re.match(r'.*([\d]:[\d]).*', karyo_string).group(1)
        fix = re.sub(r':', r';', typo)
        karyo_string = re.sub(typo, fix, karyo_string)
    karyo_string = karyo_string.replace(', ', ',')
    ## get all karyo_string before the final cell count (assuming ']' is end of karyostring)
    karyo_string = karyo_string[:karyo_string.rfind(']') + 1]
    ## start the karyostring with the pattern of <chromosome number><,>
    ## as in ....46,XX...
    m = re.match('.*?([\d]{2,3}(~[\d]{2,3})?),.*', karyo_string)
    if m:      
        karyo_string = karyo_string[m.start(1):]
    else:
        pass
        # choice here for error handling or abort for non parsable strings early on

    karyo_string = karyo_string.strip()
    return karyo_string

if __name__ == '__main__':
    print (get ('ISCN:  46,XY[20]'))
    print (get ('some 23451q23 ., example karyotype 40~43,XX,-X,add(1)(q11),der(3)del(3)(p21)t(1;3)(q21;q29),del(4)(q21),add(5)(q11.2),add(21)(q21),i(21)(q10),+1~3mar[cp4]/46,XX[16]'))