'''author@esilgard'''
#
# Copyright (c) 2015-2018 Fred Hutchinson Cancer Research Center
#
# Licensed under the Apache License, Version 2.0: http://www.apache.org/licenses/LICENSE-2.0
#
import re
import global_strings as gb
__version__ = 'iscn_parser1.0'

def get(karyotype_string):
    #print (karyotype_string)
    '''
    parse ISCN cytogenetic information from a short string of text containing
    information about the genetic variation/karyotype in the ISCN format
    ISCN formatting guidelines: http://www.cydas.org/Docs/ISCNAnalyser/Analysis.html
    list of cell types includes number of cells and dictionary of genetic abnormalities
    '''
    return_list = []
    seperate_cell_types = re.split('[/]+',karyotype_string)
    cell_type_order = 0

    for each_cell_type in seperate_cell_types:
        d = {}
        d[gb.OFFSET] = karyotype_string.find(each_cell_type) + cell_type_order
        d[gb.ABNORMALITIES] = []
        d[gb.CELL_ORDER] = cell_type_order
        cell_count = re.match(r'.*(\[c?p?([\d ]+)\]).*', each_cell_type)
        if cell_count:
            try:
                # strip off any trailing whitespace or "compositite karyotype" and cast as int
                # eventually we'll want to represent compisite counts as "no more than N"
                d[gb.CELL_COUNT] = int(cell_count.group(2).strip('cp').strip())
                each_cell_type = each_cell_type[:each_cell_type.find(cell_count.group(1))]
            ## catches error when there is no cell count
            except:
                d[gb.WARNING] = 'Error parsing number of cells'                
            # cell descr list - order is important to refer back to previous cell lines
            cell_description = each_cell_type.split(',')
            try:                
                d[gb.CHROM_NUM] = cell_description[0].strip()
                d[gb.CHROM] = cell_description[1].strip()
            except: 
                ## catches error when there is no chromosome type number and type
                d[gb.WARNING] = 'Error parsing sex chromosmes and/or number of chromosomes'
            ## if the length of the cell_description is greater than 2, there are abnormalities
            if len(cell_description) > 2:
                d[gb.ABNORMALITIES] = cell_description[2:]
                if (d[gb.CHROM] == 'sl' or 'idem' in d[gb.CHROM]) \
                and len(seperate_cell_types) > 1:
                    try:
                        d[gb.CHROM] = return_list[0][gb.CHROM]
                    except:
                        d[gb.WARNING] = 'Error parsing cell line reference to sex chromosomes'
                        d[gb.CHROM] = 'UNK'
                    d[gb.ABNORMALITIES] += return_list[0][gb.ABNORMALITIES]

                ## catch this typo - where the XX and XY is included in the idem/sl reference
                ## exclude the reference itself from the list of abnormalities
                elif (d[gb.ABNORMALITIES][0] == 'sl' or 'idem' in d[gb.ABNORMALITIES][0])\
                and len(seperate_cell_types) > 1:
                    d[gb.ABNORMALITIES] = return_list[0][gb.ABNORMALITIES] + d[gb.ABNORMALITIES][1:]
                ## when sdl is used to refer back to sl
                elif d[gb.CHROM] == 'sdl' and len(seperate_cell_types) > 2:
                    d[gb.CHROM] = return_list[1][gb.CHROM]
                    d[gb.ABNORMALITIES] += return_list[1][gb.ABNORMALITIES]
                ## catch the specific cell line references like sdl1 or sdl2 ##
                elif 'sdl' in d[gb.CHROM] or 'sl' in d[gb.CHROM]:
                    try:
                        d[gb.CHROM] = return_list[cell_type_order-1][gb.CHROM]
                        d[gb.ABNORMALITIES] = return_list[cell_type_order-1][gb.ABNORMALITIES]\
                        + d[gb.ABNORMALITIES]
                    except:      
                        d[gb.WARNING] = 'Error parsing cell line reference (e.g. sdl2)'

            ## parse abnormalities into dictionaries of type of abnormality: (num/further abn type, arm loc)
            for i in range(len(d[gb.ABNORMALITIES])):                
                if isinstance(d[gb.ABNORMALITIES][i], str):
                    # first group is the general type of abnormality 
                    # second group is the affected chromosome (sometimes involving a more complicated type of abnormality)
                    # this will get uglier eg ...   ?der(5)t(5;9)(p13;q34)t(6;9)(p22;q34)           
                    loss_gain = re.match(r'[ ]?([\+\-\?a-zA-Z ]+)[ ]?([\(]?[\d\w\~\?\(\)\-\.\;]+[\)]?)', d[gb.ABNORMALITIES][i])
                    # not doing anything with this match yet - good place to deal with polyploidy
                    copy = re.match('[xX][\d]', d[gb.ABNORMALITIES][i][-2:])                    
                    if loss_gain:   
                        # this group is the p and or q arm location (if there is one)
                        arm_location = re.match('.*([\(][\?]?[pq][pq;\?\d\.]+[\)])',loss_gain.group(0))
                        if arm_location:
                            chromosomes = loss_gain.group(2)[:loss_gain.group(2).find(arm_location.group(1))]
                            d[gb.ABNORMALITIES][i] = {loss_gain.group(1): (chromosomes, arm_location.group(1))}                            
                        else:
                            d[gb.ABNORMALITIES][i] = {loss_gain.group(1): (loss_gain.group(2), '')}     
                            
                    else:
                        ## case where there's no type of abnormality (eg del, +)
                        abnormal_chromosome = \
                        re.match(r'[ ]?([\d\-\~?\ ]*)[ ]?(mar|r|dmin|pstk+|inc)[ ]?',\
                        d[gb.ABNORMALITIES][i])
                        ## record presence of other, unrecognized abnormality
                        if abnormal_chromosome:
                            d[gb.ABNORMALITIES][i] = {'other aberration': \
                            (abnormal_chromosome.group(1), abnormal_chromosome.group(2))}
                        else:   
                            d[gb.WARNING] = 'Error finding type of abnormality'
        
        ## warning flag for error finding cell count
        else:
            d[gb.WARNING] = "Error parsing cell count"
        return_list.append(d)
        cell_type_order += 1
    return return_list

if __name__ == '__main__':
    for r in get('45,XY,-7[2]/44,idem,add(2)(q11/2),der(12)t(2;12)(q11.2;p12),-18[19]'):
        print (r)
