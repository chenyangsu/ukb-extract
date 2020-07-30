import os
import errno
import numpy as np 
import pandas as pd
import bs4 as bs
from urllib.request import urlopen
import unicodedata
from collections import OrderedDict
import argparse
# ===================================================================================================================================
# RAW DATA FILES
input_tab_file = '/path/to/ukb.tab'
hesin='/path/to/hesin.txt'
hesin_diag='/path/to/hesin_diag.txt'
coding_19_file='/path/to/coding19.tsv'
exclusion_threshold = 50 # Skip categorical feature if has more than "exclusion_threshold" categories
# ===================================================================================================================================

# CREATE OUTPUT DIRECTORY
output_dir = 'Output'
try:
    os.makedirs(output_dir)
except OSError as e:
    if e.errno != errno.EEXIST:
        raise
# ===================================================================================================================================
# SCRAPE UK BIOBANK SHOWCASE TO GET DATA-TYPE FOR EACH FIELD-ID
print('Scraping showcase webpage...')
sauce = urlopen('http://biobank.ctsu.ox.ac.uk/crystal/list.cgi').read()
soup = bs.BeautifulSoup(sauce, 'lxml')

field_listing = ['Integer',
                 'Categorical (single)',
                 'Categorical (multiple)',
                 'Continuous',
                 'Text',
                 'Date',
                 'Time',
                 'Compound']

df_descriptions = []

for url in soup.find_all('a'):
    if url.text in field_listing:
        df_type = url.text
        link = 'http://biobank.ctsu.ox.ac.uk/crystal/' +  url.get('href')
        sauce_loop = urlopen(link).read()
        soup_loop = bs.BeautifulSoup(sauce_loop, 'lxml')
        raw_text = unicodedata.normalize("NFKD", soup_loop.get_text('\t'))
        lines = raw_text.split('\n')
        for line in lines:
            splt = [s for s in line.split('\t')[1:] if s != ''
                                            and s != ' '
                                            and s!= '  '
                                            and s!= u'\u2020']
                                            #u'\u2020' is the dagger from Compound
            if len(splt) == 3 and splt[0].isdigit():
                splt.append(df_type)
                df_descriptions.append('%s\t%s\t%s\t%s'%(splt[0],
                                                         splt[1],
                                                         splt[2],
                                                         splt[3]))

with open('ukbdatafields.webscraped.tsv', 'w') as outfile:
    outfile.write('%s\n'%('\t'.join(['FieldID', 'Description', 'Category', 'ValueType'])))
    for f in df_descriptions:
        outfile.write('%s\n'%f)

# ===================================================================================================================================
# CREATE A LIST OF FIELD ID'S FOR EACH DATA-TYPE
print('Creating a list of field-ids for each data-type...')
data_types = ['Categorical_multiple', 'Categorical_single', 'Continuous', 'Integer']
for t in data_types:
    outfile = open('%s/%s.datafields'%(output_dir, t), 'w')
    with open('%s/ukbdatafields.webscraped.tsv'%output_dir, 'r') as infile:
        for line in infile:
            fid, name, category, vartype = line.rstrip().split('\t')
            # Vartype: replace brackets with an underscore
            if vartype.replace(' (', '_').replace(')','') == t:
                outfile.write('%s\n'%fid)
    outfile.close()

# ===================================================================================================================================
# SPLIT TAB_FILE BY VALUE TYPE
print('Splitting .tab file by variable type...')
fieldID_lists = [f for f in os.listdir(output_dir) if '.datafields' in f]
for fidlist in fieldID_lists:
    outprefix = fidlist.split('.')[0]
    #print('%s...'%(outprefix))

    tabfile = input_tab_file
    # GET COLUMN HEADER
    with open(tabfile, 'r') as infile:
        header = next(infile).split('\t')

    def getCols(datafields):
        cols = []
        for df in datafields:
            # Look into why you only include first visit.
            # Extract features for all visits, but make sure this doesn't break anything. 
            cols += [i for i,h in enumerate(header) if 'f.%s.0.'%df in h]
        return cols

    # UKB data-fields to search for
    with open('%s/%s'%(output_dir, fidlist), 'r') as infile:
        datafields = [line.rstrip() for line in infile]

    cols = getCols(datafields)

    with open(tabfile, 'r') as infile:
        outfile = open('%s/%s.phenoslice.csv'%(output_dir, outprefix), 'w')
        for line in infile:
            data = line.rstrip().split('\t')
            output = [data[0]] + [data[col] for col in cols]
            outfile.write('%s\n'%(','.join(output)))
        outfile.close()

# ===================================================================================================================================
# ENUMERATE POSSIBLE CATEGORIES FOR EACH FIELD ID OF TYPE 'CATEGORICAL'
print('Enumerating possible categories for each categorical-type column...')
phenoslices = ['Categorical_multiple.phenoslice.csv', 'Categorical_single.phenoslice.csv']

for phenoslice in phenoslices:
    outname = '%s.possiblecategories.txt'%phenoslice.split('.')[0]
    categories={}
    with open('%s/%s'%(output_dir, phenoslice), 'r') as infile:
        header = next(infile).rstrip().split(',')[1:]
        for f in header:
            categories[f.split('.')[1]] = set()
        for ln, line in enumerate(infile):
            #print(phenoslice, ln)
            splt, headersubset = [],[]
            for idx, i in enumerate(line.rstrip().split(',')[1:]):
                if i != 'NA':
                    splt.append(i.replace('"',''))
                    headersubset.append(header[idx])
            for i,value in enumerate(splt):
                datafield = headersubset[i].split('.')[1]
                categories[datafield].add(value)

    outfile = open('%s/%s'%(output_dir, outname), 'w')
    for i in sorted(categories):
        # Don't use spaces. There are some categories that have spaces. Ex, field 5090.
        outfile.write('%s\t%s\n'%(i, '\t'.join(sorted(categories[i]))))
    outfile.close()

# ===================================================================================================================================
# MAKE CATEGORICAL FEATURE MATRIX
print('Building categorical feature matrix...')
# NOTE: If any data-field has more than 'exclusion_threshold' number of possible categories, exclude it. 
# This will exclude data-fields like ICD9/10 codes and drug usage.

for categoricaltype in ['multiple', 'single']:
    f = 'Categorical_%s.possiblecategories.txt'%categoricaltype

    # Make a feature to index mapping 
    features = []
    with open('%s/%s'%(output_dir, f), 'r') as infile:
        for line in infile:
            splt = line.rstrip().split('\t')
            features+=['%s.%s'%(splt[0], s) for s in splt[1:]+['NA'] if len(splt) < exclusion_threshold]

    FieldDotCategoryToIndex = {x:i for i,x in enumerate(features)}
    # Define phenoslice
    phenoslice = '%s/Categorical_%s.phenoslice.csv'%(output_dir, categoricaltype)

    # Get number of individuals
    with open(phenoslice, 'r') as infile:
        n_individuals = len([i for i, line in enumerate(infile)])-1 # -1 excludes header

    # Initialize matrix
    M = np.zeros((n_individuals, len(features)))

    # Populate matrix
    with open(phenoslice, 'r') as infile:
        header = next(infile).rstrip().split(',')
        eids = []
        for rowidx,line in enumerate(infile):
            splt = line.rstrip().split(',')
            eid = splt[0]
            eids.append(eid)
            for h, category in zip(header[1:], splt[1:]):
                datafield = h.split('.')[1]
                # Create string made of current column+value seen at this individual/column
                FieldDotCategory = '%s.%s'%(datafield, category.replace('"',''))
                if FieldDotCategory in FieldDotCategoryToIndex:
                    colidx = FieldDotCategoryToIndex[FieldDotCategory]
                    M[rowidx, colidx] = 1

    # Output to file
    outname='%s.onehot.csv'%f.split('.')[0]
    with open('%s/%s'%(output_dir, outname), 'w') as outfile:
        header=','.join(['eid']+['c.%s'%f for f in features])
        outfile.write('%s\n'%header)
        for eid, row in zip(eids, M):
            outfile.write('%s,%s\n'%(eid,','.join([str(int(r)) for r in row])))

# ===================================================================================================================================
# MAKE INTEGER AND CONTINUOUS MATRICES
print('Building integer and continuous matrices...')
for phenoslice in ['%s/Integer.phenoslice.csv'%output_dir,
                   '%s/Continuous.phenoslice.csv'%output_dir]:

    prefix = os.path.basename(phenoslice).split('.')[0]

    with open(phenoslice, 'r') as infile:
        header = next(infile).rstrip().split(',')

    # Find all datafields that have more than one instance
    multiple = set()
    for h in header[1:]:
        df, nvisit, instance = h.split('.')[1:]
        if int(instance) > 0:
            multiple.add(df)

    # Get idx-list for multiple and single datafields
    single_idxs, multiple_idxs = [], []
    for i,h in enumerate(header[1:], start=1): # skip eid 
        df, nvisit, instance = h.split('.')[1:]
        if df in multiple:
            multiple_idxs.append(i)
        else:
            single_idxs.append(i)

    # Get number of individuals
    with open(phenoslice, 'r') as infile:
        n_individuals = len([i for i, line in enumerate(infile)])-1 # -1 excludes header

    # Write data to files
    f_single = open('%s/%s_single.csv'%(output_dir, prefix), 'w')
    f_multiple = open('%s/%s_multiple.csv'%(output_dir, prefix),'w')

    with open(phenoslice, 'r') as infile:
        header = next(infile).rstrip().split(',')
        f_single.write('%s\n'%(','.join([header[i] for i in [0] + single_idxs])))
        f_multiple.write('%s\n'%(','.join([header[i] for i in [0] + multiple_idxs])))

        for rowidx, line in enumerate(infile):
            data = line.rstrip().split(',')
            eid = data[0]

            if phenoslice == '%s/Integer.csv'%output_dir:
                # Did this to account for '-3: Prefer not to answer' type fields. 
                # Just set them to NA instead. This has not been tested yet. 

                dfs_single = [data[i] if '-' not in data[i] \
                              else 'NA' for i in single_idxs]

                dfs_multiple = [data[i] if '-' not in data[i] \
                                else 'NA' for i in multiple_idxs]
            else:
                dfs_single = [data[i] for i in single_idxs]
                dfs_multiple = [data[i] for i in multiple_idxs]


            f_single.write('%s,%s\n'%(eid,','.join(dfs_single)))
            f_multiple.write('%s,%s\n'%(eid,','.join(dfs_multiple)))

    # Close files
    f_single.close()
    f_multiple.close()

# ===================================================================================================================================
# EXTRACT ICD10 CODES
print('Merging hesin.txt and hesin_diag.txt tables...')

# eid and ins_index together form a primary key
hesin = pd.read_csv(hesin, sep='\t', index_col=[0,1])
# eid, ins_index and arr_index together form a primary key
hesin_diag = pd.read_csv(hesin_diag, sep='\t', index_col=[0,1,2])
# write a new table, where icd10 diagnoses have dates
df = hesin_diag.merge(hesin, how='left', on=['eid','ins_index'])
# Write merged hesin file to disk
df.to_csv('%s/hesin_merged_hesin_diag.tsv'%output_dir, sep='\t', na_rep='NA')

# In case you want to use ICD10 diagnoses that occur before baseline
baseline_visit_fieldID = 'f.53.0.0'
with open(input_tab_file, 'r') as infile:
    baseline_visit_col = next(infile).rstrip().split('\t').index(baseline_visit_fieldID)

eids = []    
baseline_visit_date = {}
with open(input_tab_file, 'r') as infile:
    next(infile)
    for line in infile:
        line = line.rstrip().split('\t', baseline_visit_col+1)
        eid = int(line[0])
        eids.append(eid)
        date = line[-2]
        baseline_visit_date[eid] = date

# Turn the multindex back into columns
df.reset_index(inplace=True) 

# Create new column 'baseline_date' in Hesin dataframe
df['baseline_visit_date'] = [baseline_visit_date[eid] for eid in df['eid']]

# Only include ICD10 codes that were diagnosed prior to baseline UKB visit. 
# To Include ICD10 codes that were diagnosed at any time, comment lines below
#df = df[df['admidate'] < df['baseline_visit_date']]
#df['baseline_date'] = [baseline_visit_date[eid] for eid in df['eid']]
#df = df[df['admidate'] < df['baseline_date']]

# --------- ICD10 LEVEL 1 --------------
print('Building ICD10 level 1 matrix...')
# Make a list of diagnoses for each individual
diagnoses = OrderedDict()
for eid in eids:
    diagnoses[eid] = []
for row in df.itertuples():
    diagnoses[row.eid].append(row.diag_icd10)

code2category = OrderedDict()
with open (coding_19_file, 'r') as infile:
    for line in infile:
        data = line.rstrip().split('\t')
        if 'Block' in data[0] or 'block' in data[0]:
            code_range = data[0].split(' ')[1]
            # Have to code these manually, bug on UKB website.
            # Ranges at ICD10 level 1 do not cover all possible codes contained at level 2. 
            # So I just manually change
            if code_range == 'M20-M25':
                code_range = 'M20-M36'
            if code_range == 'A80-A89':
                code_range = 'A80-A91'
            if code_range == 'U00-U49':
                code_range = 'U00-U81'
            if code_range == 'U82-U85':
                code_range = 'U82-U89'

            letter = code_range[0] # letter in front of digits
            lower, upper = code_range[1:3], code_range[5:7]
            codelist = ['%s%0.2d'%(letter,i) for i in range(int(lower), int(upper)+1)]
            for code in codelist:
                code2category[code] = code_range

# I do this stuff to maintain order of categories from ordered Dict. 
lookup = set()
categories = [v for k,v in code2category.items()\
              if v not in lookup and lookup.add(v) is None]

# Given a category, what is the index? 
category2index = {cat:i for i, cat in enumerate(categories)}

# output to file
outfile = open('%s/icd10level1.csv'%output_dir, 'w')
outfile.write('%s\n'%','.join(['eid']+categories))
for eid in diagnoses:
    codes = [i.replace('"', '')[:3] for i in diagnoses[eid] if i==i]
    one_hot = [0]*len(categories)
    for code in codes:
        category = code2category[code]
        index = category2index[category]
        one_hot[index] = 1
    outfile.write('%s\n'%','.join([str(eid)]+[str(o) for o in one_hot]))
outfile.close()

# --------- ICD10 LEVEL 2 --------------
print('Building ICD10 level 2 matrix...')
# Make a list of diagnoses for each individual
diagnoses = OrderedDict()
for eid in eids:
    diagnoses[eid] = []
for row in df.itertuples():
    diagnoses[row.eid].append(row.diag_icd10)

categories = []
with open(coding_19_file, 'r') as infile:
    for line in infile:
        data = line.rstrip().split('\t')
        code = data[0]
        if len(code) == 3:
            categories.append(code)

# Given a category, what is the index?
category2index = {cat:i for i, cat in enumerate(categories)}

# output to file
outfile = open('%s/icd10level2.csv'%output_dir, 'w')
outfile.write('%s\n'%','.join(['eid']+categories))

for e, eid in enumerate(diagnoses):
    #print(e)
    codes = [i.replace('"', '')[:3] for i in diagnoses[eid] if i==i]
    one_hot = [0]*len(categories)
    for code in codes:
        category = code
        index = category2index[category]
        one_hot[index] = 1
    outfile.write('%s\n'%','.join([str(eid)]+[str(o) for o in one_hot]))

# ===================================================================================================================================
# MERGE ALL FEATURES INTO A SINGLE MATRIX
print('Merging features into a single file...')
to_merge = [output_dir+'/%s'%f for f in ['Continuous_single.csv',
                                        'Continuous_multiple.csv',
                                        'Integer_single.csv',
                                        'Integer_multiple.csv',
                                        'Categorical_single.onehot.csv',
                                        'Categorical_multiple.onehot.csv', 
                                        'icd10level2.csv',
                                        'icd10level1.csv']]

outfile = open(output_dir+'/features_merged.csv', 'w')

first_file = to_merge[0]
other_files = [open(f, 'r') for f in to_merge[1:]]

with open(first_file, 'r') as infile:
    for line in infile:
        out_line = line.rstrip()
        for f in other_files:
            out_line+=(','+next(f).rstrip().split(',',1)[1])
        outfile.write('%s\n'%out_line)

for f in other_files+[outfile]:
    f.close()
