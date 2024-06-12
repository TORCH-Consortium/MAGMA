#!/usr/bin/env python3

import json
import os
import re
import argparse

import pandas as pd
import numpy as np

from tqdm import tqdm

alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
# who_v2_source = 'https://www.who.int/publications/i/item/9789240082410'

drugs = {'amikacin', 'bedaquiline', 'capreomycin', 'clofazimine', 'cycloserine', 'delamanid', 'ethambutol',
         'ethionamide', 'imipenem', 'isoniazid', 'kanamycin','levofloxacin', 'linezolid', 'meropenem', 'moxifloxacin',
         'para_aminosalicylic_acid', 'pretomanid', 'prothionamide', 'pyrazinamide', 'rifabutin', 'rifampicin',
         'streptomycin', 'terizidone'}

unknown_positions = 2.5

map_resistance_class = {
    0: 'R',
    1: 'Very high probability',
    2: 'High probability',
    3: 'Moderate probability',
    4: 'Low probability',
    5: 'Very low probability',
    6: 'S',
    -1: 'Unknown'
}

map_confidence_who = {
    'Assoc w R': 0,
    'Assoc w R - Interim': 1,
    'Uncertain significance': -1,
    'Not assoc w R - Interim': 5,
    'Not assoc w R': 6,
}

# drug_mapping = {
#     'AMI': 'amikacin',
#     'BDQ': 'bedaquiline',
#     'CAP': 'capreomycin',
#     'CFZ': 'clofazimine',
#     'DLM': 'delamanid',
#     'EMB': 'ethambutol',
#     'ETH': 'ethionamide',
#     'INH': 'isoniazid',
#     'KAN': 'kanamycin',
#     'LEV': 'levofloxacin',
#     'LZD': 'linezolid',
#     'MXF': 'moxifloxacin',
#     'PZA': 'pyrazinamide',
#     'RIF': 'rifampicin',
#     'STM': 'streptomycin',
#     'PTH': 'prothionamide'
# }
#
# drug_mapping_inv = {
#    item[1]: item[0] for item in drug_mapping.items()
# }


def create_resistance_df(sample_res, method):
    pt_df = pd.DataFrame(columns=['Drug', 'Variant', 'Resistance interpretation','Source', 'Source notation', 'Type', 'Frequency', 'Literature']).set_index(['Drug', 'Variant'])

    # DR variants from the MAGMA analysis. Add after LoFreq DR variants to give priority to MAGMA frequencies.
    for var in sample_res['dr_variants']:
        gene = var['gene_name']
        var_repr = '{}_{}'.format(gene, var['change'])

        # Add all the DR variants with classification and overwrite later if necessary
        for drug in var['gene_associated_drugs']:
            drug_name = drug.lower().replace(' ', '_')
            pt_df.loc[(drug_name, var_repr), ('Resistance interpretation', 'Source', 'Source notation', 'Type', 'Frequency', 'Literature')] = [unknown_positions, 'non-Catalogue', '', var['type'], '{:.0%}'.format(var['freq']), 'Manually curated']

        # Overwrite the variant classification for drugs which have a WHO classification last as to overrule all other classifications
        for drug in var['drugs']:
            drug_name = drug['drug'].lower().replace(' ', '_')
            pt_df.loc[(drug_name, var_repr), ('Resistance interpretation', 'Source', 'Source notation', 'Type', 'Frequency', 'Literature')] = [map_confidence_who[drug['confidence']], 'Catalogue', drug['original_mutation'], var['type'], '{:.0%}'.format(var['freq']), drug['source']]

    # Other variants from MAGMA analysis.
    for var in sample_res['other_variants']:
        gene = var['gene_name']
        var_repr = '{}_{}'.format(gene, var['change'])

        # Add all the other variants wi classification and overwrite later if necessary
        for drug in var['gene_associated_drugs']:
            drug_name = drug.lower().replace(' ', '_')
            pt_df.loc[(drug_name, var_repr), ('Resistance interpretation', 'Source', 'Source notation', 'Type', 'Frequency', 'Literature')] = [unknown_positions, 'non-Catalogue', '', var['type'], '{:.0%}'.format(var['freq']), 'Manually curated']

        # Overwrite the variant classification for drugs which have a WHO classification last as to overrule all other classifications
        if 'annotation' in var:
            for annotation in var['annotation']:
                pt_df.loc[(annotation['drug'].lower().replace(' ', '_'), var_repr), ('Resistance interpretation', 'Source', 'Source notation', 'Type', 'Frequency', 'Literature')] = [map_confidence_who[annotation['confidence']], 'Catalogue', annotation['original_mutation'], var['type'], '{:.0%}'.format(var['freq']), annotation['source']]

    pt_df['Method'] = method
    return pt_df


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Analyse resistance output from magma Pipeline')
    parser.add_argument('merged_cohort_stats_file', metavar='merged_cohort_stats_file', type=str, help='The file containing the merged cohort statistics')
    parser.add_argument('minor_res_var_dir', metavar='minor_res_var_dir', type=str, help='The directory containing the minor variants TBProfiler output files')
    parser.add_argument('struc_res_var_dir', metavar='struc_res_var_dir', type=str, help='The directory containing the structural variants TBProfiler output files')
    parser.add_argument('summary_output_dir', metavar='summary_output_dir', type=str, help='The directory where the resulting CSV files should be placed')
    args = vars(parser.parse_args())

    summary_dir = args['summary_output_dir']
    if not os.path.exists(summary_dir):
        os.makedirs(summary_dir)

    samples = {}
    for file_name in os.listdir(os.path.join(args['minor_res_var_dir'], 'results')):
        if file_name == '.DS_Store':
            continue
        keys = '.'.join(file_name.split('.')[:-2])
        samples[keys] = {}
        with open(os.path.join(args['minor_res_var_dir'], 'results', file_name)) as json_file:
            samples[keys]['lofreq'] = json.load(json_file)

    if os.path.exists(os.path.join(args['struc_res_var_dir'], 'results')):
        for file_name in os.listdir(os.path.join(args['struc_res_var_dir'], 'results')):
            if file_name == '.DS_Store':
                continue
            keys = '.'.join(file_name.split('.')[:-2])
            if keys not in samples:
                continue
            with open(os.path.join(args['struc_res_var_dir'], 'results', file_name)) as json_file:
                samples[keys]['delly'] = json.load(json_file)

#===============
# ADD FILTER FOR SAMPLES FAILING ONLY << RELABUNDANCE THRESHOLD_MET >>
#===============
    stats_df = pd.read_csv(args["merged_cohort_stats_file"], sep="\t")
    filtered_stats_df = stats_df.loc[ (stats_df["RELABUNDANCE_THRESHOLD_MET"]==0) & (stats_df["MAPPED_NTM_FRACTION_16S_THRESHOLD_MET"]==1) & (stats_df["COVERAGE_THRESHOLD_MET"]==1) & (stats_df["BREADTH_OF_COVERAGE_THRESHOLD_MET"]==1)]

    samples_df = pd.DataFrame(list(samples), columns=['full_sample'])
    filtered_samples_df = samples_df[samples_df["full_sample"].isin(filtered_stats_df["SAMPLE"].to_list())]
    samples_df = filtered_samples_df.set_index('full_sample').sort_index()


#===============
# POPULATE THE DATAFRAME
#===============
    for patient, sample in tqdm(samples_df.iterrows(), total=samples_df.shape[0]):
        sample_res = samples[patient]

        if 'lofreq' in sample_res:
            pt_df_lof = create_resistance_df(sample_res['lofreq'], 'LoFreq')
        else:
            pt_df_lof = pd.DataFrame(columns=['Drug', 'Variant']).set_index(['Drug', 'Variant'])

        if 'delly' in sample_res:
            pt_df_delly = create_resistance_df(sample_res['delly'], 'Delly')
        else:
            pt_df_delly = pd.DataFrame(columns=['Drug', 'Variant']).set_index(['Drug', 'Variant'])

        pt_df = pd.concat([pt_df_lof, pt_df_delly]).reset_index().drop_duplicates(subset=['Drug', 'Variant']).set_index(['Drug', 'Variant']).sort_index()

        for drug in list(drugs - set([i[0] for i in pt_df.index.values])):
            pt_df.loc[(drug, 'No variants found'), 'Resistance interpretation'] = 6

        pt_df['Resistance interpretation'].replace(-1, unknown_positions, inplace=True)
        for drug in pt_df.index.levels[0]:
            conc = min(pt_df.loc[drug, 'Resistance interpretation'].value_counts().keys())
            pt_df.loc[drug, 'Conclusion'] = conc

        pt_df = pt_df.reset_index().sort_values(['Conclusion', 'Drug', 'Resistance interpretation', 'Variant'])
        for column in ['Conclusion', 'Resistance interpretation']:
            pt_df[column] = pt_df[column].apply(lambda c: map_resistance_class[-1] if c == unknown_positions else None if pd.isna(c) else map_resistance_class[c])

        pt_df = pt_df[['Drug', 'Conclusion', 'Variant', 'Resistance interpretation', 'Type', 'Frequency', 'Method', 'Literature', 'Source notation']]

        # Write the sheet to excel with formatting
        with pd.ExcelWriter(os.path.join(summary_dir, '{}.xlsx'.format(patient)), engine='xlsxwriter') as writer:
            pt_df.set_index(['Drug', 'Conclusion', 'Variant']).to_excel(writer, sheet_name='resistance_variants')  # Updated sheet name

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
                writer.sheets['resistance_variants'].conditional_format('{}1:{}{}'.format(alphabet[1], alphabet[1], pt_df.shape[0]+1),  cond_format)
                writer.sheets['resistance_variants'].conditional_format('{}1:{}{}'.format(alphabet[3], alphabet[3], pt_df.shape[0]+1),  cond_format)

            # Autofit the worksheet and hide columns
            # writer.sheets['resistance_variants'].autofit()  # Not available in the version we are using
            writer.sheets['resistance_variants'].set_column(2, 2, 25, None, None)
            writer.sheets['resistance_variants'].set_column(5, 5, 25, None, None)
