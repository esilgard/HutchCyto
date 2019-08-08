'''author@esilgard'''
#
# Copyright (c) 2015-2018 Fred Hutchinson Cancer Research Center
#
# Licensed under the Apache License, Version 2.0: http://www.apache.org/licenses/LICENSE-2.0
#

import global_strings as gb
import eln_classification
import re

__version__='cytogenetics_abnormality_parser1.0'

def get(cell_list, karyotype_string, clonality_flag):
    '''
        take a parsed list of cells, where each element within the cell group 
        has a list of dictionaries of abnormalities as well as a cell count, 
        a chromosome number, a warning flag, a chromosome, and a cell order 
    '''
    if not clonality_flag :
        min_abn = 0; min_mono = 0
    else:
        min_abn = 2; min_mono = 3
    return_dictionary_list = [{gb.NAME:gb.KARYO_STR, gb.VALUE:karyotype_string, 
        gb.VERSION:__version__}]
   
    # chromosomes that can be generally coded with the six coding groups:
    # translocation(chr,arm), -chr, add(chr,arm), del(chr,arm), inv(chr,arm),dup or trp(chr,arm) 
    all_chromosomes = [str(i) for i in range(1,22)] + ['X','Y']
    other_chr_group = ['13','14','15','21','22']  # no short arm abnormalities
    arm_list = ['p','q']

    # specific translocations of interest in AML risk classification
    specific_translocations = ['15;17', '6;9', '14;16', '9;22','8;21', '16;16',
        '4;14', '11;14', '16;20', '3;3']
    ## a dictionary of abnormality types and their cell counts
    abnormalities = {}
    offsets = {}
    
    def add_to_d(abnormality, cell_count, variation_start, variation_end):
        '''
        helper function to put the abnormalities cell counts and the 
        start and stop positions of the abnormality string into the 
        abnormalities and offsets dictionaries
        '''
        if cell_count:
            # None cell_counts are binary values; only integers get added
            abnormalities[abnormality] = abnormalities.get(abnormality, 0) + cell_count
        elif abnormality == gb.ABNS:
            abnormalities[abnormality] = 0
        else:
            abnormalities[abnormality] = True
        offsets[abnormality] = offsets.get(abnormality, [])
        offsets[abnormality].append((variation_start, variation_end)) 

    all_monosomy_set = set([])
    all_trisomy_set = set([])
    all_abnormality_set = set([])
    ###########################################################################
    abnormalities[gb.WARNING] = False; offsets[gb.WARNING] = []
    abnormalities[gb.MONO_TYPE] = False; offsets[gb.MONO_TYPE] = []
    abnormalities[gb.CMPX_TYPE] = False; offsets[gb.CMPX_TYPE] = []
    ## start by counting cells with each type of pertinent aberration      
    for x in cell_list:    
        if x.get(gb.WARNING, False):
            abnormalities[gb.WARNING] = 'Error carried over from iscn_parser'            
        try:
            # used to track unique abnormalities
            abnormality_set = set([])
            monosomy_set = set([])
            other_structural_abnormalities_set = set([])
            cell_count = x[gb.CELL_COUNT]           
            cell_offset = x[gb.OFFSET]
            ## NOTE - currently this only gets the lower limit of cells in cases of "43-47"
            chromosome_number = int(x[gb.CHROM_NUM][:2])
            if x[gb.ABNORMALITIES]: 
                for y in x[gb.ABNORMALITIES]:
                    try:  
                        for z, zz in y.items():       
                            variation_string = z + zz[0] + zz[1]
                            stripped_chr = zz[0].strip('(').strip(')')
                            if cell_count >= min_abn:
                                variation_start = cell_offset + \
                                (karyotype_string[cell_offset:].find(variation_string))
                                # note - this does NOT capture character offsets correctly
                                # for 'idem', 'sl'  references
                                variation_end = variation_start + len(variation_string)
                                #....should this just be the all_abnormality set?
                                #add_to_d(gb.ABNS, 1, variation_start, variation_end)
                                ## all trisomies
                                if z in ['+','-']:
                                    all_abnormality_set.add(variation_string) 
                                    abnormality_set.add(variation_string)
                                    if z == '+':
                                        try:
                                            # only add in "full" trisomies
                                            int(stripped_chr)
                                            all_trisomy_set.add(variation_string)
                                            add_to_d(gb.TRIS, 1, variation_start, variation_end)
                                        except:
                                            pass                                    
                                        if stripped_chr in all_chromosomes:  
                                            add_to_d('+' + stripped_chr, cell_count, variation_start, variation_end)
                                    ## all monosomies                            
                                    elif z == '-' and cell_count >= min_mono: 
                                        # track only autosomal monosomies for 
                                        # monosomal karyotype classification
                                        if stripped_chr not in ['X','Y']:
                                            monosomy_set.add(variation_string)
                                            all_monosomy_set.add(variation_string)
                                            add_to_d(gb.MONOS, 1, variation_start, variation_end)
                                        if stripped_chr in all_chromosomes:   
                                            add_to_d('-' + stripped_chr, cell_count, variation_start, variation_end)
                                    # track other structural abnormalities for 
                                    # monosomal karyotype classificaiton
                                else:
                                    if z not in ['r','mar','+ma','+mar','add']:
                                        other_structural_abnormalities_set.add(variation_string) 
                                    all_abnormality_set.add(variation_string) 
                                    abnormality_set.add(variation_string)
                                # specific salient translocations
                                if z == 't' or 'dic' in z or z == '?t':
                                    if stripped_chr in specific_translocations:
                                        add_to_d('t(' + stripped_chr + ')', cell_count, variation_start, variation_end)\
                                    # general translocations involving p or q arms for encoding group
                                    # just q arm for the "other" chromosome group
                                    for arm in arm_list:
                                        for each in all_chromosomes:
                                            if each in other_chr_group and arm == 'p':
                                                pass
                                            else:
                                                chr_list = stripped_chr.split(';')
                                                if each in chr_list:
                                                    location = chr_list.index(each)
                                                    if arm in zz[1].split(';')[location]:
                                                        add_to_d('translocation(' + each + arm + ')', cell_count, variation_start, variation_end)
                                   
                                # derivative chromosomes (generally more complicated strings)
                                elif 'der' in z:
                                    for arm in arm_list:
                                        for each in all_chromosomes:
                                            if each in other_chr_group and arm == 'p':
                                                pass
                                            else:
                                                chr_list = stripped_chr.split(';')                                        
                                                # if der involves other abnormality, won't necessarily 
                                                #split chromosomes cleanly e.g.  ['12)t(12','15']
                                                for chr_part in chr_list:                                                    
                                                    if each in chr_part:                                                        
                                                        location = chr_list.index(chr_part)                                                       
                                                        # a derivitive from q10 on will mean a full loss of the p arm
                                                        other_arm = [a for a in arm_list if a!=arm][0]
                                                        if arm + '10' in zz[1].split(';')[location]:
                                                            add_to_d('del(' + each + other_arm + ')', cell_count, variation_start, variation_end)
                                                        else:
                                                            if arm in  zz[1].split(';')[location]:
                                                                add_to_d('translocation(' + each + arm + ')', cell_count, variation_start, variation_end)

                                ## isochromes - implicit deletion of p or q arms
                                elif z == 'i' or z == '?i':
                                    if stripped_chr in all_chromosomes:
                                        if 'p' in zz[1]:
                                            add_to_d('del(' + stripped_chr + 'q)', cell_count, variation_start, variation_end)                                            
                                        elif 'q' in zz[1]:  
                                            add_to_d('del(' + stripped_chr + 'p)', cell_count, variation_start, variation_end)
                                    
                                else:
                                    ## re.search allows for variants with '?'
                                    generic_abn_label = re.match('^.{0,5}(dup|trp|inv|ins|del|add).{0,5}$',z)
                                    if generic_abn_label:
                    
                                        label = generic_abn_label.group(1)                                        
                                        for arm in arm_list:
                                            if stripped_chr in all_chromosomes and arm in zz[1]:    
                                                if stripped_chr in other_chr_group and arm == 'p':                                                 
                                                    pass
                                                else:                                        
                                                    if label == 'ins': label = 'translocation'                                                    
                                                    add_to_d(label + '('+ stripped_chr + arm + ')', cell_count, variation_start, variation_end)
                                    
                    ## catch any other formatting abnormalities/parsing errors
                    except:
                        abnormalities[gb.WARNING] = True; x[gb.WARNING] = 'unknown parsing err'
                                                
            ## there are no abnormalities - add up the "normal" cells
            elif x[gb.CHROM] in ['XX','XY']:
                add_to_d(gb.NORMAL, cell_count, cell_offset, cell_offset + 2)
            else:
                ## other sec chromosome abnormalities
                add_to_d(gb.SEX_CHRM_ABN, cell_count, cell_offset, cell_offset + len(x[gb.CHROM]))
            ## complex and monosomal karyotype classifications (no char offsets currently)            
            if len(abnormality_set) > 2:
                abnormalities[gb.CMPX_TYPE] = True
            if len(monosomy_set) > 1 or (len(monosomy_set) == 1 and \
                len(other_structural_abnormalities_set)) > 0:
                abnormalities[gb.MONO_TYPE] = True            
        ## catch trouble with cell counts etc       
        except:
            abnormalities[gb.WARNING] = True; x[gb.WARNING] = 'err w/ cell counts'   
    add_to_d(gb.ABNS, len(all_abnormality_set), 0, len(karyotype_string))                  
    #append all abnormality info to return_list
    for each_variation in abnormalities:
        offsets[each_variation] = offsets.get(each_variation, [])
        return_dictionary_list.append({gb.NAME:each_variation, \
            gb.VALUE:abnormalities[each_variation],gb.VERSION:__version__})
   
    return_dictionary_list.append(eln_classification.get(abnormalities, abnormality_set, offsets, karyotype_string, 0))

    return return_dictionary_list