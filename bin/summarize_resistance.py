#!/usr/bin/env python3

import json
import os
import re
import argparse

import pandas as pd

from tqdm import tqdm

alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
drugs = set(['amikacin', 'bedaquiline', 'capreomycin', 'clofazimine', 'cycloserine', 'delamanid', 'ethambutol', 'ethionamide', 'imipenem', 'isoniazid', 'isoniazid_high_dose', 'kanamycin','levofloxacin', 'linezolid', 'meropenem', 'moxifloxacin', 'moxifloxacin_high_dose', 'para_aminosalicylic_acid', 'pretomanid', 'prothionamide', 'pyrazinamide', 'rifabutin', 'rifampicin', 'rifampicin_high_dose', 'streptomycin', 'terizidone'])
#drugs = set(['amikacin', 'bedaquiline', 'capreomycin', 'clofazimine', 'delamanid', 'ethambutol', 'ethionamide', 'isoniazid', 'isoniazid_high_dose', 'kanamycin', 'levofloxacin', 'linezolid', 'moxifloxacin', 'prothionamide', 'pyrazinamide', 'rifampicin', 'streptomycin'])
class_map = {
    0: 'R',
    1: 'U',
    2: 'S'
}

major_variants_sheet_name = 'Major variants'
minor_variants_sheet_name = 'Minor variants'

def extract_patient_and_sample(full_sample_name):
    if 'SMARTT' in full_sample_name:
        return pd.Series(['-'.join(full_sample_name.split('-')[:3]), '-'.join(full_sample_name.split('-')[3:])])
    if 'FSDOH' in full_sample_name:
        return pd.Series(['-'.join(full_sample_name.split('-')[:2]), '-'.join(full_sample_name.split('-')[2:])])
    return pd.Series([full_sample_name, full_sample_name])

def add_var_to_df(df, pt_id, drug, var, freq):
    if drug not in list(df) or pd.isna(df.loc[patient, drug]):
        df.loc[pt_id, drug] = '{} ({:.0%})'.format(var_repr, freq)
    else:
        if var_repr not in df.loc[pt_id, drug]:
            df.loc[pt_id, drug] += ' & {} ({:.0%})'.format(var_repr, freq)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Analyse resistance output from XBS Pipeline')
    parser.add_argument('major_res_var_dir', metavar='major_res_var_dir', type=str, help='The direcotry containing the major variants TBProfiler output files')
    parser.add_argument('minor_res_var_dir', metavar='minor_res_var_dir', type=str, help='The direcotry containing the minor variants TBProfiler output files')
    parser.add_argument('summary_output_dir', metavar='summary_output_dir', type=str, help='The directory where the resulting excel sheets should be placed')
    args = vars(parser.parse_args())

    summary_dir = args['summary_output_dir']
    if not os.path.exists(summary_dir):
        os.makedirs(summary_dir)

    samples = {}
    for file_name in os.listdir(os.path.join(args['major_res_var_dir'], 'results')):
        keys = file_name.split('.')
        samples[keys[1]] = {}
        with open(os.path.join(os.path.join(args['major_res_var_dir'], 'results', file_name))) as json_file:
            samples[keys[1]]['xbs'] = json.load(json_file)

    if os.path.exists(os.path.join(os.path.join(args['minor_res_var_dir'], 'results'))):
        for file_name in os.listdir(os.path.join(os.path.join(args['minor_res_var_dir'], 'results'))):
            keys = file_name.split('.')
            if keys[1] not in samples:
                continue
                #samples[keys[1]] = {}
            with open(os.path.join(os.path.join(os.path.join(args['minor_res_var_dir'], 'results', file_name)))) as json_file:
                samples[keys[1]]['lofreq'] = json.load(json_file)

    samples_df = pd.DataFrame(list(samples), columns=['full_sample'])
    #samples_df[['patient', 'sample']] = samples_df['full_sample'].apply(lambda sample: extract_patient_and_sample(sample))
    samples_df = samples_df.set_index('full_sample').sort_index()

    for patient, sample in tqdm(samples_df.iterrows(), total=samples_df.shape[0]):
        sample_res = samples[patient]

        pt_df_xbs = pd.DataFrame(columns=['Drug', 'Variant', 'Interpretation', 'Source'] + ['Conclusion {}'.format(patient)] + list([patient])).set_index(['Drug', 'Variant'])
        pt_df_lof = pd.DataFrame(columns=['Drug', 'Variant', 'Interpretation', 'Source'] + ['Conclusion {}'.format(patient)] + list([patient])).set_index(['Drug', 'Variant'])

        """
        Add the DR variants from the XBS analysis to the xbs variant dataframe.
        Do this after adding the lofreq dr variants to show the XBS variant frequencies.
        """
        for var in sample_res['xbs']['dr_variants']:
            gene = var['gene']
            if gene == '.':
                gene = var['locus_tag']
            var_repr = '{}_{}'.format(gene, var['change'])
            for drug in var['drugs']:
                drug = drug['drug'].lower().replace(' ', '_')
                pt_df_xbs.loc[(drug, var_repr), (patient, 'Interpretation', 'Source')] = ['{:.0%}'.format(var['freq']), 0, 'WHO Catalogue']

        """
        Add the other variants from the XBS analysis to the xbs variant dataframe.
        """
        for var in sample_res['xbs']['other_variants']:
            gene = var['gene']
            if gene == '.':
                gene = var['locus_tag']
            var_repr = '{}_{}'.format(gene, var['change'])
            # Add all the other variants as unknown classification and overwrite their classification later if necessary
            for drug in var['gene_associated_drugs']:
                drug = drug.lower().replace(' ', '_')
                pt_df_xbs.loc[(drug, var_repr), (patient, 'Interpretation', 'Source')] = ['{:.0%}'.format(var['freq']), 1, 'Tier 1 or 2 gene']
            # Overwrite the variant classification for drugs which have a WHO sens classification last as to overrule all other classifications
            if 'annotation' in var:
                for annotation in var['annotation']:
                    if annotation['type'] == 'resistance_association_confidence' and (int(annotation['confidence']) == 4 or int(annotation['confidence']) == 5):
                        pt_df_xbs.loc[(annotation['drug'].lower().replace(' ', '_'), var_repr), (patient, 'Interpretation', 'Source')] = ['{:.0%}'.format(var['freq']), 2, 'WHO Catalogue']
                    elif annotation['type'] == 'resistance_association_confidence' and (int(annotation['confidence']) == 3):
                        pt_df_xbs.loc[(annotation['drug'].lower().replace(' ', '_'), var_repr), (patient, 'Interpretation', 'Source')] = ['{:.0%}'.format(var['freq']), 1, 'WHO Catalogue']
                    else:
                        display(annotation)

        """
        Add the DR variants from the lofreq analysis to the xbs variant dataframe.
        """
        if 'lofreq' in sample_res:
            for var in sample_res['lofreq']['dr_variants']:
                gene = var['gene']
                if gene == '.':
                    gene = var['locus_tag']
                var_repr = '{}_{}'.format(gene, var['change'])
                for drug in var['drugs']:
                    drug = drug['drug'].lower().replace(' ', '_')
                    pt_df_lof.loc[(drug, var_repr), (patient, 'Interpretation', 'Source')] = ['{:.0%}'.format(var['freq']), 0, 'WHO Catalogue']

        """
        Add the other variants from the lofreq analysis to the lofreq variant dataframe.
        """
        if 'lofreq' in sample_res:
            for var in sample_res['lofreq']['other_variants']:
                gene = var['gene']
                if gene == '.':
                    gene = var['locus_tag']
                var_repr = '{}_{}'.format(gene, var['change'])
                # Add all the other variants as unknown classification and overwrite their classification later if necessary
                for drug in var['gene_associated_drugs']:
                    drug = drug.lower().replace(' ', '_')
                    pt_df_lof.loc[(drug, var_repr), (patient, 'Interpretation', 'Source')] = ['{:.0%}'.format(var['freq']), 1, 'Tier 1 or 2 gene']
                # Overwrite the variant classification for drugs which have a WHO sens classification last as to overrule all other classifications
                if 'annotation' in var:
                    for annotation in var['annotation']:
                        if annotation['type'] == 'resistance_association_confidence' and (int(annotation['confidence']) == 4 or int(annotation['confidence']) == 5):
                            pt_df_lof.loc[(annotation['drug'].lower().replace(' ', '_'), var_repr), (patient, 'Interpretation', 'Source')] = ['{:.0%}'.format(var['freq']), 2, 'WHO Catalogue']
                        elif annotation['type'] == 'resistance_association_confidence' and (int(annotation['confidence']) == 3):
                            pt_df_lof.loc[(annotation['drug'].lower().replace(' ', '_'), var_repr), (patient, 'Interpretation', 'Source')] = ['{:.0%}'.format(var['freq']), 1, 'WHO Catalogue']
                        else:
                            display(annotation)

        # Remove all variants in the XBS summary from the lofreq summary
        pt_df_lof = pt_df_lof.drop([i for i in pt_df_xbs.index if i in pt_df_lof.index])
        for drug in list(drugs - set([i[0] for i in pt_df_xbs.index.values])):
            pt_df_xbs.loc[(drug, 'No variants found'), ('Interpretation')] = 2
        for drug in list(drugs - set([i[0] for i in pt_df_lof.index.values])):
            pt_df_lof.loc[(drug, 'No variants found'), ('Interpretation')] = 2
        
        
        """
        Calculate the conclusion for all samples
        """
        for _, sample in samples_df.loc[[patient]].iterrows():
            work_df = pt_df_xbs[pd.notna(pt_df_xbs[patient])]
            for drug in pt_df_xbs.index.levels[0]:
                if drug not in set([i[0] for i in work_df.index.values]):
                    pt_df_xbs.loc[drug, 'Conclusion {}'.format(patient)] = 2
                elif 0 in work_df.loc[drug, 'Interpretation'].value_counts():
                    pt_df_xbs.loc[drug, 'Conclusion {}'.format(patient)] = 0
                elif 1 in work_df.loc[drug, 'Interpretation'].value_counts():
                    pt_df_xbs.loc[drug, 'Conclusion {}'.format(patient)] = 1
                else:
                    pt_df_xbs.loc[drug, 'Conclusion {}'.format(patient)] = 2
        for _, sample in samples_df.loc[[patient]].iterrows():
            work_df = pt_df_lof[pd.notna(pt_df_lof[patient])]
            for drug in pt_df_lof.index.levels[0]:
                if drug not in set([i[0] for i in work_df.index.values]):
                    pt_df_lof.loc[drug, 'Conclusion {}'.format(patient)] = 2
                elif 0 in work_df.loc[drug, 'Interpretation'].value_counts():
                    pt_df_lof.loc[drug, 'Conclusion {}'.format(patient)] = 0
                elif 1 in work_df.loc[drug, 'Interpretation'].value_counts():
                    pt_df_lof.loc[drug, 'Conclusion {}'.format(patient)] = 1
                else:
                    pt_df_lof.loc[drug, 'Conclusion {}'.format(patient)] = 2

        x2 = pt_df_lof.copy(deep=True)
        pt_df_xbs = pt_df_xbs.reset_index().sort_values(['Conclusion {}'.format(patient)] + ['Drug', 'Interpretation', 'Variant'])
        for column in [i for i in pt_df_xbs if 'Conclusion' in i or 'Interpretation' == i]:
            pt_df_xbs[column] = pt_df_xbs[column].apply(lambda c: class_map[c])
        pt_df_lof = pt_df_lof.reset_index().sort_values(['Conclusion {}'.format(patient)] + ['Drug', 'Interpretation', 'Variant'])
        for column in [i for i in pt_df_lof if 'Conclusion' in i or 'Interpretation' == i]:
            pt_df_lof[column] = pt_df_lof[column].apply(lambda c: class_map[c])

        """
        Write both sheets to excel with formatting"""
        with pd.ExcelWriter(os.path.join(summary_dir, '{}.xlsx'.format(patient)), engine='xlsxwriter') as writer:
            pt_df_xbs.set_index(['Drug'] + ['Conclusion {}'.format(patient)] + ['Variant']).to_excel(writer, sheet_name=major_variants_sheet_name)
            pt_df_lof.set_index(['Drug'] + ['Conclusion {}'.format(patient)] + ['Variant']).to_excel(writer, sheet_name=minor_variants_sheet_name)
            #unclassified.reset_index().sort_values(by=['Conclusion', 'Drug', 'Interpretation', 'Variant']).set_index(['Drug', 'Conclusion', 'Variant']).to_excel(writer, sheet_name='Unclassified Variants')

            # Create some formatting
            format_sens = writer.book.add_format({'bold': False, 'font_color': 'green'})
            format_res = writer.book.add_format({'bold': True, 'font_color': 'red'})
            cond_res = {'type': 'cell', 'criteria': '==', 'value': '"R"', 'format': format_res}
            cond_sens = {'type': 'cell', 'criteria': '==', 'value': '"S"', 'format': format_sens}

            # Add formatting to Variants sheet
            for i in range(samples_df.loc[[patient]].shape[0]):
                writer.sheets[major_variants_sheet_name].conditional_format('{}1:{}{}'.format(alphabet[1+i], alphabet[1+i], pt_df_xbs.shape[0]+1),  cond_res)
                writer.sheets[major_variants_sheet_name].conditional_format('{}1:{}{}'.format(alphabet[1+i], alphabet[1+i], pt_df_xbs.shape[0]+1),  cond_sens)
            writer.sheets[major_variants_sheet_name].conditional_format('{}1:{}{}'.format(alphabet[2+samples_df.loc[[patient]].shape[0]], alphabet[2+samples_df.loc[[patient]].shape[0]], pt_df_xbs.shape[0]+1),  cond_res)
            writer.sheets[major_variants_sheet_name].conditional_format('{}1:{}{}'.format(alphabet[2+samples_df.loc[[patient]].shape[0]], alphabet[2+samples_df.loc[[patient]].shape[0]], pt_df_xbs.shape[0]+1),  cond_sens)

            # Add formatting to Lofreq variants sheet
            for i in range(samples_df.loc[[patient]].shape[0]):
                writer.sheets[minor_variants_sheet_name].conditional_format('{}1:{}{}'.format(alphabet[1+i], alphabet[1+i], pt_df_lof.shape[0]+1),  cond_res)
                writer.sheets[minor_variants_sheet_name].conditional_format('{}1:{}{}'.format(alphabet[1+i], alphabet[1+i], pt_df_lof.shape[0]+1),  cond_sens)
            writer.sheets[minor_variants_sheet_name].conditional_format('{}1:{}{}'.format(alphabet[2+samples_df.loc[[patient]].shape[0]], alphabet[2+samples_df.loc[[patient]].shape[0]], pt_df_lof.shape[0]+1),  cond_res)
            writer.sheets[minor_variants_sheet_name].conditional_format('{}1:{}{}'.format(alphabet[2+samples_df.loc[[patient]].shape[0]], alphabet[2+samples_df.loc[[patient]].shape[0]], pt_df_lof.shape[0]+1),  cond_sens)
