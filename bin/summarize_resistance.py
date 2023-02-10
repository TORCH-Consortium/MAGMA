#!/usr/bin/env python3

import json
import os
import re
import argparse

import warnings
warnings.filterwarnings("ignore")

import pandas as pd

from tqdm import tqdm

alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
drugs = set(['amikacin', 'bedaquiline', 'capreomycin','clofazimine', 'cycloserine', 'delamanid', 'ethambutol', 'ethionamide', 'imipenem', 'isoniazid', 'isoniazid_high_dose', 'kanamycin','levofloxacin', 'linezolid', 'meropenem', 'moxifloxacin', 'moxifloxacin_high_dose', 'para_aminosalicylic_acid', 'pretomanid', 'prothionamide', 'pyrazinamide', 'rifabutin', 'rifampicin', 'rifampicin_high_dose', 'streptomycin', 'terizidone'])
class_map = {
    0: 'R',
    1: 'U',
    2: 'S'
}

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
    parser.add_argument('xbs_output_dir', metavar='XBS_output_dir', type=str, help='The input sample file (the study id cannot start with \'XBS_REF_\')')
    parser.add_argument('xbs_run_id', metavar='XBS_run_id', type=str, help='The directory to which all output files should be written')
    parser.add_argument('literature_db', metavar='literature_db', type=str, help='The directory to which all output files should be written')
    parser.add_argument('literature_db_forced', metavar='literature_db_forced', type=str, help='The directory to which all output files should be written')
    parser.add_argument('summary_output_dir', metavar='summary_output_dir', type=str, help='The name of the output VCF file')

    args = vars(parser.parse_args())

    run = args['xbs_run_id']
    xbs_output_dir = args['xbs_output_dir']
    summary_dir = os.path.join(args['summary_output_dir'], run)
    if not os.path.exists(summary_dir):
        os.makedirs(summary_dir)
    lit_df = pd.read_csv(args['literature_db'], index_col=[0,1])
    lit_df_forced = pd.read_csv(args['literature_db_forced'], index_col=[0,1])

    samples = {}
    for file_name in os.listdir(os.path.join(xbs_output_dir, run, 'analyses/drug_resistance/major_variants/results')):
    #for file_name in os.listdir(os.path.join(xbs_output_dir, 'resistance', run,'XBS/results')):
        keys = file_name.split('.')
#        if 'SMARTT' not in keys[0]:
#            continue
        samples[keys[1]] = {}
        with open(os.path.join(os.path.join(xbs_output_dir, run, 'analyses/drug_resistance/major_variants/results', file_name))) as json_file:
            samples[keys[1]]['xbs'] = json.load(json_file)

    if os.path.exists(os.path.join(os.path.join(xbs_output_dir, run, 'analyses/drug_resistance/minor_variants/results'))):
        for file_name in os.listdir(os.path.join(xbs_output_dir, run, 'analyses/drug_resistance/minor_variants/results')):
        #for file_name in os.listdir(os.path.join(xbs_output_dir, 'resistance', run,'lofreq/results')):
            keys = file_name.split('.')
#            if 'SMARTT' not in keys[0]:
#                continue
            if keys[1] not in samples:
                continue
                samples[keys[1]] = {}
            with open(os.path.join(os.path.join(xbs_output_dir, run, 'analyses/drug_resistance/minor_variants/results', file_name))) as json_file:
                samples[keys[1]]['lofreq'] = json.load(json_file)

    samples_df = pd.DataFrame(list(samples), columns=['full_sample'])
    samples_df[['patient', 'sample']] = samples_df['full_sample'].apply(lambda sample: extract_patient_and_sample(sample))
    samples_df = samples_df.set_index('patient').sort_index()

    for patient in tqdm(sorted(set(samples_df.index))):
        pt_df_xbs = pd.DataFrame(columns=['Drug', 'Variant', 'Interpretation', 'Source'] + ['Conclusion {}'.format(i) for i in list(samples_df.loc[[patient]]['full_sample'])] + list(samples_df.loc[[patient]]['full_sample'])).set_index(['Drug', 'Variant'])
        pt_df_lof = pd.DataFrame(columns=['Drug', 'Variant', 'Interpretation', 'Source'] + ['Conclusion {}'.format(i) for i in list(samples_df.loc[[patient]]['full_sample'])] + list(samples_df.loc[[patient]]['full_sample'])).set_index(['Drug', 'Variant'])

        """
        Generate the resistance summary files for all samples
        """
        for _, sample in samples_df.loc[[patient]].iterrows():
            sample_res = samples[sample['full_sample']]

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
                        pt_df_xbs.loc[(drug, var_repr), (sample['full_sample'], 'Interpretation', 'Source')] = ['{:.0%}'.format(var['freq']), 0, 'WHO Catalogue']
                    if var_repr in lit_df_forced.index:
                        for drug, classification in lit_df_forced.loc[var_repr]['classification'].iteritems():
                            pt_df_xbs.loc[(drug.lower().replace(' ', '_'), var_repr), (sample['full_sample'], 'Interpretation', 'Source')] = ['{:.0%}'.format(var['freq']), classification, 'Experts']

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
                    pt_df_xbs.loc[(drug, var_repr), (sample['full_sample'], 'Interpretation', 'Source')] = ['{:.0%}'.format(var['freq']), 0, 'WHO Catalogue']
                if var_repr in lit_df_forced.index:
                    for drug, classification in lit_df_forced.loc[var_repr]['classification'].iteritems():
                        pt_df_xbs.loc[(drug.lower().replace(' ', '_'), var_repr), (sample['full_sample'], 'Interpretation', 'Source')] = ['{:.0%}'.format(var['freq']), classification, 'Experts']

            """
            Add the other variants from the lofreq analysis to the xbs variant dataframe.
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
                        pt_df_lof.loc[(drug, var_repr), (sample['full_sample'], 'Interpretation', 'Source')] = ['{:.0%}'.format(var['freq']), 1, '']
                    # Overwrite the variant classification for drugs which have a WHO sens classification last as to overrule all other classifications
                    if 'annotation' in var:
                        for annotation in var['annotation']:
                            if annotation['type'] == 'resistance_association_confidence' and (int(annotation['confidence']) == 4 or int(annotation['confidence']) == 5):
                                pt_df_lof.loc[(annotation['drug'].lower().replace(' ', '_'), var_repr), (sample['full_sample'], 'Interpretation', 'Source')] = ['{:.0%}'.format(var['freq']), 2, 'WHO Catalogue']
                            elif annotation['type'] == 'resistance_association_confidence' and (int(annotation['confidence']) == 3):
                                pt_df_lof.loc[(annotation['drug'].lower().replace(' ', '_'), var_repr), (sample['full_sample'], 'Interpretation', 'Source')] = ['{:.0%}'.format(var['freq']), 1, 'WHO Catalogue']
                            else:
                                display(annotation)
                    # Overwrite the variant classification if it was previously classified by the experts
                    if var_repr in lit_df.index:
                        for drug, classification in lit_df.loc[var_repr]['classification'].iteritems():
                            pt_df_lof.loc[(drug.lower().replace(' ', '_'), var_repr), (sample['full_sample'], 'Interpretation', 'Source')] = ['{:.0%}'.format(var['freq']), classification, 'Experts']
                    if var_repr in lit_df_forced.index:
                        for drug, classification in lit_df_forced.loc[var_repr]['classification'].iteritems():
                            pt_df_lof.loc[(drug.lower().replace(' ', '_'), var_repr), (sample['full_sample'], 'Interpretation', 'Source')] = ['{:.0%}'.format(var['freq']), classification, 'Experts']

            """
            Add the other variants from the XBS analysis to the xbs variant dataframe.
            Do this after adding the lofreq other variants to show the XBS variant frequencies.
            """
            for var in sample_res['xbs']['other_variants']:
                gene = var['gene']
                if gene == '.':
                    gene = var['locus_tag']
                var_repr = '{}_{}'.format(gene, var['change'])
                # Add all the other variants as unknown classification and overwrite their classification later if necessary
                for drug in var['gene_associated_drugs']:
                    drug = drug.lower().replace(' ', '_')
                    pt_df_xbs.loc[(drug, var_repr), (sample['full_sample'], 'Interpretation', 'Source')] = ['{:.0%}'.format(var['freq']), 1, '']
                # Overwrite the variant classification for drugs which have a WHO sens classification last as to overrule all other classifications
                if 'annotation' in var:
                    for annotation in var['annotation']:
                        if annotation['type'] == 'resistance_association_confidence' and (int(annotation['confidence']) == 4 or int(annotation['confidence']) == 5):
                            pt_df_xbs.loc[(annotation['drug'].lower().replace(' ', '_'), var_repr), (sample['full_sample'], 'Interpretation', 'Source')] = ['{:.0%}'.format(var['freq']), 2, 'WHO Catalogue']
                        elif annotation['type'] == 'resistance_association_confidence' and (int(annotation['confidence']) == 3):
                            pt_df_xbs.loc[(annotation['drug'].lower().replace(' ', '_'), var_repr), (sample['full_sample'], 'Interpretation', 'Source')] = ['{:.0%}'.format(var['freq']), 1, 'WHO Catalogue']
                        else:
                            display(annotation)
                # Overwrite the variant classification if it was previously classified by the experts
                if var_repr in lit_df.index:
                    for drug, classification in lit_df.loc[var_repr]['classification'].iteritems():
                        pt_df_xbs.loc[(drug.lower().replace(' ', '_'), var_repr), (sample['full_sample'], 'Interpretation', 'Source')] = ['{:.0%}'.format(var['freq']), classification, 'Experts']
                if var_repr in lit_df_forced.index:
                    for drug, classification in lit_df_forced.loc[var_repr]['classification'].iteritems():
                        pt_df_xbs.loc[(drug.lower().replace(' ', '_'), var_repr), (sample['full_sample'], 'Interpretation', 'Source')] = ['{:.0%}'.format(var['freq']), classification, 'Experts']

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
            work_df = pt_df_xbs[pd.notna(pt_df_xbs[sample['full_sample']])]
            for drug in pt_df_xbs.index.levels[0]:
                if drug not in set([i[0] for i in work_df.index.values]):
                    pt_df_xbs.loc[drug, 'Conclusion {}'.format(sample['full_sample'])] = 2
                elif 0 in work_df.loc[drug, 'Interpretation'].value_counts():
                    pt_df_xbs.loc[drug, 'Conclusion {}'.format(sample['full_sample'])] = 0
                elif 1 in work_df.loc[drug, 'Interpretation'].value_counts():
                    pt_df_xbs.loc[drug, 'Conclusion {}'.format(sample['full_sample'])] = 1
                else:
                    pt_df_xbs.loc[drug, 'Conclusion {}'.format(sample['full_sample'])] = 2
        for _, sample in samples_df.loc[[patient]].iterrows():
            work_df = pt_df_lof[pd.notna(pt_df_lof[sample['full_sample']])]
            for drug in pt_df_lof.index.levels[0]:
                if drug not in set([i[0] for i in work_df.index.values]):
                    pt_df_lof.loc[drug, 'Conclusion {}'.format(sample['full_sample'])] = 2
                elif 0 in work_df.loc[drug, 'Interpretation'].value_counts():
                    pt_df_lof.loc[drug, 'Conclusion {}'.format(sample['full_sample'])] = 0
                elif 1 in work_df.loc[drug, 'Interpretation'].value_counts():
                    pt_df_lof.loc[drug, 'Conclusion {}'.format(sample['full_sample'])] = 1
                else:
                    pt_df_lof.loc[drug, 'Conclusion {}'.format(sample['full_sample'])] = 2

        x2 = pt_df_lof.copy(deep=True)
        pt_df_xbs = pt_df_xbs.reset_index().sort_values(['Conclusion {}'.format(i) for i in samples_df.loc[[patient], 'full_sample']] + ['Drug', 'Interpretation', 'Variant'])
        for column in [i for i in pt_df_xbs if 'Conclusion' in i or 'Interpretation' == i]:
            pt_df_xbs[column] = pt_df_xbs[column].apply(lambda c: class_map[c])
        pt_df_lof = pt_df_lof.reset_index().sort_values(['Conclusion {}'.format(i) for i in samples_df.loc[[patient], 'full_sample']] + ['Drug', 'Interpretation', 'Variant'])
        for column in [i for i in pt_df_lof if 'Conclusion' in i or 'Interpretation' == i]:
            pt_df_lof[column] = pt_df_lof[column].apply(lambda c: class_map[c])

        """
        Write both sheets to excel with formatting"""
        with pd.ExcelWriter(os.path.join(summary_dir, '{}.xlsx'.format(patient)), engine='xlsxwriter') as writer:
            pt_df_xbs.set_index(['Drug'] + ['Conclusion {}'.format(i) for i in samples_df.loc[[patient], 'full_sample']] + ['Variant']).to_excel(writer, sheet_name='Variants')
            pt_df_lof.set_index(['Drug'] + ['Conclusion {}'.format(i) for i in samples_df.loc[[patient], 'full_sample']] + ['Variant']).to_excel(writer, sheet_name='Lofreq variants')
            #unclassified.reset_index().sort_values(by=['Conclusion', 'Drug', 'Interpretation', 'Variant']).set_index(['Drug', 'Conclusion', 'Variant']).to_excel(writer, sheet_name='Unclassified Variants')

            # Create some formatting
            format_sens = writer.book.add_format({'bold': False, 'font_color': 'green'})
            format_res = writer.book.add_format({'bold': True, 'font_color': 'red'})
            cond_res = {'type': 'cell', 'criteria': '==', 'value': '"R"', 'format': format_res}
            cond_sens = {'type': 'cell', 'criteria': '==', 'value': '"S"', 'format': format_sens}

            # Add formatting to Variants sheet
            for i in range(samples_df.loc[[patient]].shape[0]):
                writer.sheets['Variants'].conditional_format('{}1:{}{}'.format(alphabet[1+i], alphabet[1+i], pt_df_xbs.shape[0]+1),  cond_res)
                writer.sheets['Variants'].conditional_format('{}1:{}{}'.format(alphabet[1+i], alphabet[1+i], pt_df_xbs.shape[0]+1),  cond_sens)
            writer.sheets['Variants'].conditional_format('{}1:{}{}'.format(alphabet[2+samples_df.loc[[patient]].shape[0]], alphabet[2+samples_df.loc[[patient]].shape[0]], pt_df_xbs.shape[0]+1),  cond_res)
            writer.sheets['Variants'].conditional_format('{}1:{}{}'.format(alphabet[2+samples_df.loc[[patient]].shape[0]], alphabet[2+samples_df.loc[[patient]].shape[0]], pt_df_xbs.shape[0]+1),  cond_sens)

            # Add formatting to Lofreq variants sheet
            for i in range(samples_df.loc[[patient]].shape[0]):
                writer.sheets['Lofreq variants'].conditional_format('{}1:{}{}'.format(alphabet[1+i], alphabet[1+i], pt_df_lof.shape[0]+1),  cond_res)
                writer.sheets['Lofreq variants'].conditional_format('{}1:{}{}'.format(alphabet[1+i], alphabet[1+i], pt_df_lof.shape[0]+1),  cond_sens)
            writer.sheets['Lofreq variants'].conditional_format('{}1:{}{}'.format(alphabet[2+samples_df.loc[[patient]].shape[0]], alphabet[2+samples_df.loc[[patient]].shape[0]], pt_df_lof.shape[0]+1),  cond_res)
            writer.sheets['Lofreq variants'].conditional_format('{}1:{}{}'.format(alphabet[2+samples_df.loc[[patient]].shape[0]], alphabet[2+samples_df.loc[[patient]].shape[0]], pt_df_lof.shape[0]+1),  cond_sens)
