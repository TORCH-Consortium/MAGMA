#!/usr/bin/env python3

import json
import os
import re
import argparse

import pandas as pd
import numpy as np

from tqdm import tqdm

ALPHABET = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
DRUGS = {'amikacin', 'bedaquiline', 'capreomycin', 'clofazimine', 'cycloserine', 'delamanid', 'ethambutol',
         'ethionamide', 'imipenem', 'isoniazid', 'kanamycin', 'levofloxacin', 'linezolid', 'meropenem', 'moxifloxacin',
         'para_aminosalicylic_acid', 'pretomanid', 'prothionamide', 'pyrazinamide', 'rifabutin', 'rifampicin',
         'streptomycin', 'terizidone'}
UNKNOWN_POSITIONS = 2.5

MAJOR_VARIANTS_SHEET_NAME = 'resistance_variants'

MAP_RESISTANCE_CLASS = {
    0: 'R',
    1: 'Very high probability',
    2: 'High probability',
    3: 'Moderate Probability',
    4: 'Low probability',
    5: 'Very low probability',
    6: 'S',
    -1: 'Unknown'
}

MAP_CONFIDENCE_WHO = {
    1: 0,
    2: 1,
    3: -1,
    4: 5,
    5: 6,
}

MAP_CONFIDENCE_WHO_NAMING = {
    1: 'Associated with resistance',
    2: 'Associated with resistance - interim',
    3: 'Uncertain significance',
    4: 'Not associated with resistance - interim',
    5: 'Not associated with resistance',
}

DRUG_MAPPING = {
    'AMI': 'amikacin',
    'BDQ': 'bedaquiline',
    'CAP': 'capreomycin',
    'CFZ': 'clofazimine',
    'DLM': 'delamanid',
    'EMB': 'ethambutol',
    'ETH': 'ethionamide',
    'INH': 'isoniazid',
    'KAN': 'kanamycin',
    'LEV': 'levofloxacin',
    'LZD': 'linezolid',
    'MXF': 'moxifloxacin',
    'PZA': 'pyrazinamide',
    'RIF': 'rifampicin',
    'STM': 'streptomycin',
    'PTH': 'prothionamide'
}

drug_mapping_inv = {
   item[1]: item[0] for item in DRUG_MAPPING.items()
}

def create_resistance_df(sample_res, method='XBS'):
    pt_df = pd.DataFrame(columns=['Drug', 'Variant', 'Resistance interpretation','WHO Catalogue', 'Source', 'Source notation', 'Type', 'Frequency', 'Literature', 'Observations', 'Fraction']).set_index(['Drug', 'Variant'])
    """
    Add the DR variants from the magma analysis to the magma variant dataframe.
    Do this after adding the lofreq dr variants to show the magma variant frequencies.
    """
    for var in sample_res['dr_variants']:
        gene = var['gene_name']
        if gene == '.':
            gene = var['locus_tag']
        var_repr = '{}_{}'.format(gene, var['change'])

        #Add all drugs for which this variant falls in a tier 1 or 2 gene
        for drug in var['gene_associated_drugs']:
            drug_name = drug.lower().replace(' ', '_')
            pt_df.loc[(drug, var_repr), ('Frequency', 'Resistance interpretation', 'WHO Catalogue', 'Type', 'Source')] = ['{:.0%}'.format(var['freq']), -1, 'Not listed', var['type'], 'Tier 1 or 2 gene']

        # Overwrite the unknown for the resistant drugs
        for drug in var['drugs']:
            drug_name = drug['drug'].lower().replace(' ', '_')
            pt_df.loc[(drug_name, var_repr), ('Frequency', 'Resistance interpretation', 'WHO Catalogue', 'Type', 'Source', 'Source notation', 'Literature')] = ['{:.0%}'.format(var['freq']), MAP_CONFIDENCE_WHO[int(drug['confidence'])], MAP_CONFIDENCE_WHO_NAMING[int(drug['confidence'])], var['type'], 'Catalogue', drug['who original'], drug['literature']]

        # Overwrite the unknown for the sensitive drugs
        for drug in var['annotation']:
            drug_name = drug['drug'].lower().replace(' ', '_')
            pt_df.loc[(drug_name, var_repr), ('Frequency', 'Resistance interpretation', 'WHO Catalogue', 'Type', 'Source', 'Source notation', 'Literature')] = ['{:.0%}'.format(var['freq']), MAP_CONFIDENCE_WHO[int(drug['who_confidence'])], MAP_CONFIDENCE_WHO_NAMING[int(drug['who_confidence'])], var['type'], 'Catalogue', drug['who_original'], drug['literature']]
            if float(drug['who_solo_res']) + float(drug['who_sens']) == 0:
                pt_df.loc[(drug_name, var_repr), ('Observations')] = ['{}R - {}S'.format(int(float(drug['who_solo_res'])), int(float(drug['who_sens'])))]
            else:
                pt_df.loc[(drug_name, var_repr), ('Observations', 'Fraction')] = ['{}R - {}S'.format(int(float(drug['who_solo_res'])), int(float(drug['who_sens']))), '{:.2f}'.format(float(drug['who_r_fraction']))]
                    

    """
    Add the other variants from the magma analysis to the magma variant dataframe.
    """
    for var in sample_res['other_variants']:
        gene = var['gene_name']
        if gene == '.':
            gene = var['locus_tag']
        var_repr = '{}_{}'.format(gene, var['change'])

        # Add all the other variants as unknown classification and overwrite their classification later if necessary
        for drug in var['gene_associated_drugs']:
            drug = drug.lower().replace(' ', '_')
            pt_df.loc[(drug, var_repr), ('Frequency', 'Resistance interpretation', 'WHO Catalogue', 'Type', 'Source')] = ['{:.0%}'.format(var['freq']), -1, 'Not listed', var['type'], 'Tier 1 or 2 gene']

        # Overwrite the variant classification for drugs which have a WHO sens classification last as to overrule all other classifications
        if 'annotation' in var:
            for annotation in var['annotation']:
                if annotation['type'] == 'who_confidence':
                    if int(annotation['who_confidence']) != 3:
                        pt_df.loc[(annotation['drug'].lower().replace(' ', '_'), var_repr), ('Observations', 'Fraction')] = None
                    if int(annotation['who_confidence']) > 3:
                        pt_df.loc[(annotation['drug'].lower().replace(' ', '_'), var_repr), ('Frequency', 'Resistance interpretation', 'WHO Catalogue', 'Type', 'Source', 'Source notation', 'Literature')] = ['{:.0%}'.format(var['freq']), MAP_CONFIDENCE_WHO[int(annotation['who_confidence'])], MAP_CONFIDENCE_WHO_NAMING[int(annotation['who_confidence'])], var['type'], 'Catalogue', annotation['who_original'], annotation['literature']]
                    elif int(annotation['who_confidence']) == 3:
                        if float(annotation['who_solo_res']) + float(annotation['who_sens']) == 0:
                            pt_df.loc[(annotation['drug'].lower().replace(' ', '_'), var_repr), ('Frequency', 'Resistance interpretation', 'WHO Catalogue', 'Type', 'Source', 'Source notation', 'Literature', 'Observations')] = ['{:.0%}'.format(var['freq']), MAP_CONFIDENCE_WHO[int(annotation['who_confidence'])], MAP_CONFIDENCE_WHO_NAMING[int(annotation['who_confidence'])], var['type'], 'Catalogue', annotation['who_original'], annotation['literature'], '{}R - {}S'.format(int(float(annotation['who_solo_res'])), int(float(annotation['who_sens'])))]
                        else:
                            pt_df.loc[(annotation['drug'].lower().replace(' ', '_'), var_repr), ('Frequency', 'Resistance interpretation', 'WHO Catalogue', 'Type', 'Source', 'Source notation', 'Literature', 'Observations', 'Fraction')] = ['{:.0%}'.format(var['freq']), MAP_CONFIDENCE_WHO[int(annotation['who_confidence'])], MAP_CONFIDENCE_WHO_NAMING[int(annotation['who_confidence'])], var['type'], 'Catalogue', annotation['who_original'], annotation['literature'], '{}R - {}S'.format(int(float(annotation['who_solo_res'])), int(float(annotation['who_sens']))), '{:.2f}'.format(float(annotation['who_r_fraction']))]
                    else:
                        display(annotation)
                else:
                    display(annotation)
    pt_df['Method'] = method
    return pt_df

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Analyse resistance output from magma Pipeline')
    parser.add_argument('merged_cohort_stats_file', metavar='merged_cohort_stats_file', type=str, help='The file containing the merged cohort statistics')
    parser.add_argument('major_res_var_dir', metavar='major_res_var_dir', type=str, help='The directory containing the major variants TBProfiler output files')
    parser.add_argument('minor_res_var_dir', metavar='minor_res_var_dir', type=str, help='The directory containing the minor variants TBProfiler output files')
    parser.add_argument('struc_res_var_dir', metavar='struc_res_var_dir', type=str, help='The directory containing the structural variants TBProfiler output files')
    parser.add_argument('summary_output_dir', metavar='summary_output_dir', type=str, help='The directory where the resulting excel sheets should be placed')
    args = vars(parser.parse_args())


    summary_dir = args['summary_output_dir']
    if not os.path.exists(summary_dir):
        os.makedirs(summary_dir)

    samples = {}
    for file_name in os.listdir(os.path.join(args['major_res_var_dir'], 'results')):
        if file_name == '.DS_Store':
            continue
        keys = '.'.join(file_name.split('.')[:-2])
        samples[keys] = {}
        with open(os.path.join(os.path.join(args['major_res_var_dir'], 'results', file_name))) as json_file:
            samples[keys]['magma'] = json.load(json_file)

    if os.path.exists(os.path.join(os.path.join(args['minor_res_var_dir'], 'results'))):
        for file_name in os.listdir(os.path.join(os.path.join(args['minor_res_var_dir'], 'results'))):
            if file_name == '.DS_Store':
                continue
            keys = '.'.join(file_name.split('.')[:-2])
            if keys not in samples:
                continue
                #samples[keys[1]] = {}
            with open(os.path.join(os.path.join(os.path.join(args['minor_res_var_dir'], 'results', file_name)))) as json_file:
                samples[keys]['lofreq'] = json.load(json_file)

    if os.path.exists(os.path.join(os.path.join(args['struc_res_var_dir'], 'results'))):
        for file_name in os.listdir(os.path.join(os.path.join(args['struc_res_var_dir'], 'results'))):
            if file_name == '.DS_Store':
                continue
            keys = '.'.join(file_name.split('.')[:-2])
            if keys not in samples:
                continue
                #samples[keys[1]] = {}
            with open(os.path.join(os.path.join(os.path.join(args['struc_res_var_dir'], 'results', file_name)))) as json_file:
                samples[keys]['delly'] = json.load(json_file)

#===============
# ADD FILTER FOR SAMPLES PASSING << ALL_THRESHOLDS >>
#===============
    stats_df = pd.read_csv(args["merged_cohort_stats_file"], sep="\t")
    filtered_stats_df = stats_df[ stats_df["ALL_THRESHOLDS_MET"] == 1 ]

    samples_df = pd.DataFrame(list(samples), columns=['full_sample'])
    #samples_df[['patient', 'sample']] = samples_df['full_sample'].apply(lambda sample: extract_patient_and_sample(sample))
    filtered_samples_df = samples_df[samples_df["full_sample"].isin(filtered_stats_df["SAMPLE"].to_list())]
    samples_df = filtered_samples_df.set_index('full_sample').sort_index()

#===============
# POPULATE THE DATAFRAME
#===============

    for patient, sample in tqdm(samples_df.iterrows(), total=samples_df.shape[0]):
        sample_res = samples[patient]

        pt_df_magma = create_resistance_df(sample_res['magma'], 'XBS')

        if 'lofreq' in sample_res:
            pt_df_lof = create_resistance_df(sample_res['lofreq'], 'LoFreq')
        else:
            pt_df_lof = pd.DataFrame(columns=['Drug', 'Variant']).set_index(['Drug', 'Variant'])

        if 'delly' in sample_res:
            pt_df_delly = create_resistance_df(sample_res['delly'], 'Delly')
        else:
            pt_df_delly = pd.DataFrame(columns=['Drug', 'Variant']).set_index(['Drug', 'Variant'])

        pt_df = pd.concat([pt_df_magma, pt_df_lof, pt_df_delly]).reset_index().drop_duplicates(subset=['Drug', 'Variant']).set_index(['Drug', 'Variant']).sort_index()

        for drug in list(DRUGS - set([i[0] for i in pt_df.index.values])):
            pt_df.loc[(drug, 'No variants found'), ('WHO Catalogue')] = ''
            pt_df.loc[(drug, 'No variants found'), ('Resistance interpretation')] = 6

        pt_df['Resistance interpretation'].replace(-1, UNKNOWN_POSITIONS, inplace=True)
        for drug in pt_df.index.levels[0]:
            conc = min(pt_df.loc[drug, 'Resistance interpretation'].value_counts().keys())
            pt_df.loc[drug, 'Conclusion'] = conc

        pt_df = pt_df.reset_index().sort_values(['Conclusion', 'Drug', 'Resistance interpretation', 'Variant'])
        for column in ['Conclusion', 'Resistance interpretation']:
            pt_df[column] = pt_df[column].apply(lambda c: MAP_RESISTANCE_CLASS[-1] if c == UNKNOWN_POSITIONS else None if pd.isna(c) else MAP_RESISTANCE_CLASS[c])

        pt_df = pt_df[['Drug', 'Conclusion', 'Variant', 'Resistance interpretation', 'WHO Catalogue', 'Type', 'Frequency', 'Method', 'Literature', 'Source notation', 'Observations', 'Fraction']]

        """
        Write the sheet to excel with formatting
        """
        with pd.ExcelWriter(os.path.join(summary_dir, '{}.xlsx'.format(patient)), engine='xlsxwriter') as writer:
            pt_df.set_index(['Drug', 'Conclusion', 'Variant']).to_excel(writer, sheet_name=MAJOR_VARIANTS_SHEET_NAME)

            format_sens = writer.book.add_format({'bold': False, 'font_color': 'green'})
            format_mod = writer.book.add_format({'bold': True, 'font_color': 'orange'})
            format_res = writer.book.add_format({'bold': True, 'font_color': 'red'})
            cond_res_0 = {'type': 'cell', 'criteria': '==', 'value': '"R"', 'format': format_res}
            cond_res_1 = {'type': 'cell', 'criteria': '==', 'value': '"Very high probability"', 'format': format_res}
            cond_res_2 = {'type': 'cell', 'criteria': '==', 'value': '"High probability"', 'format': format_res}
            cond_mod_3 = {'type': 'cell', 'criteria': '==', 'value': '"Moderate probability"', 'format': format_mod}
            cond_sens_4 = {'type': 'cell', 'criteria': '==', 'value': '"Low probability"', 'format': format_sens}
            cond_sens_5 = {'type': 'cell', 'criteria': '==', 'value': '"Very low probability"', 'format': format_sens}
            cond_sens_6 = {'type': 'cell', 'criteria': '==', 'value': '"S"', 'format': format_sens}

            # Add formatting to Variants sheet
            for cond_format in [cond_res_0, cond_res_1, cond_res_2, cond_mod_3, cond_sens_4, cond_sens_5, cond_sens_6]:
                writer.sheets[MAJOR_VARIANTS_SHEET_NAME].conditional_format('{}1:{}{}'.format(ALPHABET[1], ALPHABET[1], pt_df.shape[0] + 1), cond_format)
                writer.sheets[MAJOR_VARIANTS_SHEET_NAME].conditional_format('{}1:{}{}'.format(ALPHABET[3], ALPHABET[3], pt_df.shape[0] + 1), cond_format)

            # Autofit the worksheet and hide columns
            writer.sheets[MAJOR_VARIANTS_SHEET_NAME].autofit() # Not available in the version we are using
            writer.sheets[MAJOR_VARIANTS_SHEET_NAME].set_column(2, 2, 25, None, None)
            writer.sheets[MAJOR_VARIANTS_SHEET_NAME].set_column(5, 5, 25, None, None)
            #writer.sheets[major_variants_sheet_name].set_column(5, 6, None, None, {"hidden": True})
