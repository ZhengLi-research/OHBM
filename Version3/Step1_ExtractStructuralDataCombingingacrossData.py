# The function of this script is to extract structural data from freesurfer
################ add global structural index!!!!!!!!!!!!!!!!!!!!!
# prepare data for plotting and GAM in R
# 为每个数据集的每个结构指标生成一个tsv文件。没行代表每个人，每列代表一个结构指标、以及人口统计学指标

import numpy as np
import glob
import pandas as pd

datapath = '/Users/lizheng/Desktop/RBC/'

qc_version =  'pass' # 'artifact' or 'qc'

# Generate tsv files for RBC structural data

dtype = ['ct', 'sa', 'mc', 'gc', 'fi', 'ci', 'gv', 'sv'] # ct: cortex thickness; sa: surface area; mc: mean curvature; gc: gaussian curvature; fi: folding index; ci: curvature index; gv: gray matter volume
dataset = ['hbn'] #'bhrc', 'ccnp', 'nki', 'pnc'

for iDataset in dataset: # 按照数据库类型进行循环
    #if qc_version == 'noqc':
    #    fileNames = sorted(glob.glob(datapath + 'data/%s_Freesurfer/'
    #              '*_regionsurfacestats.tsv'
    #              % iDataset.upper()))
    # elif qc_version == 'artifact':
    #     fileNames = sorted(glob.glob(datapath + 'data/%s_Freesurfer_artifact/'
    #                 '*_regionsurfacestats.tsv'
    #                 % iDataset.upper()))
    if qc_version == 'pass':
        fileNames = sorted(glob.glob((datapath + '%s/%s_Freesurfer/freesurfer/sub-*/'
                                      '*_regionsurfacestats.tsv') 
                                      % (iDataset.upper(), iDataset.upper())))
        
                                                      
    subj_list = [iFile.split('/')[-1].split('_')[0].split('-')[1]
                 for iFile in fileNames]
    
    # extract specific structural features from tabulated freesurfer data
    # cortex thickness
    ct = []
    # surface area
    sa = []
    # mean curvature
    mc = []
    # gaussian curvature
    gc = []
    # folding index
    fi = []
    # curvature index
    ci = []
    # gray matter volume
    gv = []
    # subcortical structures
    sv = []


    for i, iFile in enumerate(fileNames):
        temp_file = iFile
        temp_file2 = iFile
        brain_data = pd.read_csv(iFile, delimiter='\t')
        DK68 = brain_data.loc[brain_data['atlas'] == 'glasser', ('hemisphere', 'StructName', 'SurfArea', 'GrayVol', 'ThickAvg', 'MeanCurv', 'GausCurv', 'FoldInd', 'CurvInd')]
        DK68.reset_index(inplace=True, drop=True)
        temp_filefin = temp_file2.replace('regionsurfacestats', 'brainmeasures')
        brain_global_file = pd.read_csv(temp_filefin, delimiter='\t')
        brain_global_index = brain_global_file[['EstimatedTotalIntraCranialVol_eTIV']] #'Cortex_PialSurfArea_lh', 'Cortex_PialSurfArea_rh', 'Cortex_MeanThickness_lh', 'Cortex_MeanThickness_rh'
        subcortical_file = brain_global_file.copy()
        #sv_index = subcortical_file[['Left_Thalamus_Proper_Volume_mm3', 'Right_Thalamus_Proper_Volume_mm3', 'Left_Caudate_Volume_mm3', 'Right_Caudate_Volume_mm3', 'Left_Putamen_Volume_mm3', 'Right_Putamen_Volume_mm3'
        #                    , 'Left_Pallidum_Volume_mm3', 'Right_Pallidum_Volume_mm3', 'Left_Hippocampus_Volume_mm3', 'Right_Hippocampus_Volume_mm3', 'Left_Amygdala_Volume_mm3', 'Right_Amygdala_Volume_mm3'
        #                    , 'Left_Accumbens_area_Volume_mm3', 'Right_Accumbens_area_Volume_mm3']]
        # 合并字符串
        DK68['StructName'] = DK68['hemisphere'] + '_' + DK68['StructName'] 
        region_names = list(DK68['StructName'])
        # sv_names = list(sv_index.columns)
        region_names.extend(['EstimatedTotalIntraCranialVol_eTIV']) # 'Cortex_PialSurfArea_lh', 'Cortex_PialSurfArea_rh', 'Cortex_MeanThickness_lh', 'Cortex_MeanThickness_rh', 
        # sv_names.extend(['EstimatedTotalIntraCranialVol_eTIV'])
    
        ct.append(np.concatenate([np.array(DK68['ThickAvg']), np.array(brain_global_index).flatten()]))# ct.append(np.array(DK68['ThickAvg']))
        sa.append(np.concatenate([np.array(DK68['SurfArea']), np.array(brain_global_index).flatten()]))# sa.append(np.array(DK68['SurfArea']))
        mc.append(np.concatenate([np.array(DK68['MeanCurv']), np.array(brain_global_index).flatten()]))# mc.append(np.array(DK68['MeanCurv']))
        gc.append(np.concatenate([np.array(DK68['GausCurv']), np.array(brain_global_index).flatten()]))# gc.append(np.array(DK68['GausCurv']))
        fi.append(np.concatenate([np.array(DK68['FoldInd']), np.array(brain_global_index).flatten()]))# fi.append(np.array(DK68['FoldInd']))
        ci.append(np.concatenate([np.array(DK68['CurvInd']), np.array(brain_global_index).flatten()]))# ci.append(np.array(DK68['CurvInd']))
        gv.append(np.concatenate([np.array(DK68['GrayVol']), np.array(brain_global_index).flatten()]))# gv.append(np.array(DK68['GrayVol']))
        # sv.append(np.concatenate([np.array(sv_index).flatten(), np.array(brain_global_index['EstimatedTotalIntraCranialVol_eTIV'])]))# sv.append(np.array(sv_index).flatten())
        print('\nFile %i/%i done!' % (i, len(fileNames)))

    for iType in dtype:
        # load demographics and qc data
        demogs = pd.read_csv((datapath + '%s/%s_BIDS/' +
                      'study-%s_desc-participants.tsv') 
                     % (iDataset.upper(), iDataset.upper(), iDataset.upper()),
                     delimiter='\t')
        demogs['participant_id'] = np.array(demogs['participant_id']).astype(
            'str')
        
        qc_data = pd.read_csv((datapath + '%s/%s_FreeSurfer/' +
                            'study-%s_desc-T1_qc.tsv') 
                            % (iDataset.upper(), iDataset.upper(), iDataset.upper()),
                            delimiter='\t')
        qc_data['participant_id'] = np.array(qc_data['participant_id']).astype(
            'str')
        
        # generate metric dataframe
        if iType == 'ct':
            df_metric = pd.DataFrame(np.array(ct), columns=region_names)
            df_metric['participant_id'] = subj_list
            df_metric['meanVal'] = np.mean(np.array(ct), axis=1)
        elif iType == 'gv':
            df_metric = pd.DataFrame(np.array(gv), columns=region_names)
            df_metric['participant_id'] = subj_list
            df_metric['meanVal'] = np.mean(np.array(gv), axis=1)
        elif iType == 'sa':
            df_metric = pd.DataFrame(np.array(sa), columns=region_names)
            df_metric['participant_id'] = subj_list
            df_metric['meanVal'] = np.mean(np.array(sa), axis=1)
        elif iType == 'mc':
            df_metric = pd.DataFrame(np.array(mc), columns=region_names)
            df_metric['participant_id'] = subj_list
            df_metric['meanVal'] = np.mean(np.array(mc), axis=1)
        elif iType == 'gc':
            df_metric = pd.DataFrame(np.array(gc), columns=region_names)
            df_metric['participant_id'] = subj_list
            df_metric['meanVal'] = np.mean(np.array(gc), axis=1)
        elif iType == 'fi':
            df_metric = pd.DataFrame(np.array(fi), columns=region_names)
            df_metric['participant_id'] = subj_list
            df_metric['meanVal'] = np.mean(np.array(fi), axis=1)
        elif iType == 'ci':
            df_metric = pd.DataFrame(np.array(ci), columns=region_names)
            df_metric['participant_id'] = subj_list
            df_metric['meanVal'] = np.mean(np.array(ci), axis=1)
        #elif iType == 'sv':
        #    df_metric = pd.DataFrame(np.array(sv), columns=sv_names)
        #    df_metric['participant_id'] = subj_list
        #    df_metric['meanVal'] = np.mean(np.array(ci), axis=1)
        # if data is 'bhrc' or 'nki', only keep baseline scans
        if iDataset == 'bhrc':
            demogs_ses1 = demogs.query('session_id == 1')
            demogs_ses1.reset_index(inplace=True, drop=True)
            del demogs
            demogs = demogs_ses1.copy()
            del demogs_ses1

            qc_data_ses1 = qc_data.query('session_id == 1')
            qc_data_ses1.reset_index(inplace=True, drop=True)
            del qc_data
            qc_data = qc_data_ses1.copy()
            del qc_data_ses1

        elif iDataset == 'nki':
            demogs_ses1 = demogs[demogs['session_id'].str.contains('BAS1')]
            demogs_ses1.reset_index(inplace=True, drop=True)
            del demogs
            demogs = demogs_ses1.copy()
            del demogs_ses1

            qc_data_ses1 = qc_data[qc_data['session_id'].str.contains('BAS1')]
            qc_data_ses1.reset_index(inplace=True, drop=True)
            del qc_data
            qc_data = qc_data_ses1.copy()
            del qc_data_ses1      
        # add demographics
        x = np.array(df_metric['participant_id'])
        y = np.array(demogs['participant_id'])
        xy, x_ind, y_ind = np.intersect1d(x, y, return_indices=True)

        df_demogs = df_metric.iloc[x_ind, :]
        df_demogs.reset_index(inplace=True, drop=True)

        demogs_shared = demogs.iloc[y_ind, :]
        demogs_shared.reset_index(inplace=True, drop=True)

        df_demogs_shared = pd.merge(df_demogs, demogs_shared,
                                    on='participant_id', how='left')

        # add T1 qc info
        x = np.array(df_demogs_shared['participant_id'])
        y = np.array(qc_data['participant_id'])
        xy, x_ind, y_ind = np.intersect1d(x, y, return_indices=True)

        df_qc = df_demogs_shared.iloc[x_ind, :]
        df_qc.reset_index(inplace=True, drop=True)

        qc_data_shared = qc_data.iloc[y_ind, :]
        qc_data_shared.reset_index(inplace=True, drop=True)

        df_qc['euler'] = qc_data_shared['euler'].values
        df_qc['qc_determination'] = qc_data_shared['qc_determination'].values

        # check whether everyone has euler and age values
        # check_nan_euler = df_qc['euler'].isnull().values.any()
        # check_nan_age = df_qc['age'].isnull().values.any()
        # if check_nan_euler or check_nan_age:
        #     euler_nanidx = np.where(np.isnan(df_qc['euler'].values))
        #     age_nanidx = np.where(np.isnan(df_qc['age'].values))
        #     all_nanidx = list(age_nanidx[0]) + list(euler_nanidx[0])
        #     df_qc.drop(np.array(all_nanidx), inplace=True)
        #     df_qc.reset_index(inplace=True, drop=True)                 
        
        # extract data for age range 6-22 years old
        df_final = df_qc.copy()
        df_final = df_final.query('5 <= age <= 22')
        df_final.reset_index(inplace=True, drop=True)
        df_final.drop(columns=['lh_???', 'rh_???'], inplace=True)
        


        # if qc_version == 'noqc':
        #     output_filename = (datapath + 'data/dataR/%s_df_%s_noqc.tsv'
        #                        % (iDataset, iType))
        # elif qc_version == 'artifact':
        #     output_filename = (datapath + 'data/dataR/%s_df_%s_artifact.tsv'
        #                        % (iDataset, iType))
        if qc_version == 'pass':
            output_filename = ( '/Users/lizheng/Desktop/同步文件夹/博士研究课题/OHBM会议数据分析/Version3/Structural_data/%s_df_%s_pass.tsv'
                               % (iDataset, iType))           

        df_final.to_csv(output_filename, sep='\t')

# combine all 5 datasets

# dtype = ['ct_pass', 'gv_pass', 'sa_pass', 'mc_pass', 'gc_pass', 'fi_pass', 'ci_pass']

# for iType in dtype:
#     pnc = pd.read_csv(datapath + 'data/dataR/pnc_df_%s.tsv' % iType,
#                       delimiter='\t')
#     hbn = pd.read_csv(datapath + 'data/dataR/hbn_df_%s.tsv' % iType,
#                       delimiter='\t')
#     bhrc = pd.read_csv(datapath + 'data/dataR/bhrc_df_%s.tsv' % iType,
#                        delimiter='\t')
#     nki = pd.read_csv(datapath + 'data/dataR/nki_df_%s.tsv' % iType,
#                       delimiter='\t')
#     ccnp = pd.read_csv(datapath + 'data/dataR/ccnp_df_%s.tsv' % iType,
#                        delimiter='\t')

#     combined_df = pd.concat([pnc, hbn, bhrc, nki, ccnp], axis=0)

#     combined_df.reset_index(inplace=True, drop=True)

#     combined_df.to_csv(datapath + 'data/dataR/combined_df_%s.tsv' % iType,
#                        index=False, sep='\t')