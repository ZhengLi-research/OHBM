%%% Script that creates input files to plsgui
%%% IV 2023
clear, clc
datapath='/Users/lizheng/Desktop/同步文件夹/博士研究课题/OHBM会议数据分析/Version3/Structural_data/Structural_split/';
behaviorpath = '/Users/lizheng/Desktop/同步文件夹/博士研究课题/OHBM会议数据分析/Version3/output/';
plspath='/Users/lizheng/Desktop/同步文件夹/博士研究课题/OHBM会议数据分析/Version3/output/'; %% test path output
% /Users/lizheng/Desktop/同步文件夹/博士研究课题/OHBM会议数据分析/Version2/output/

%% Prep datafile
% data is indexed by datafile, mat-file that reads into matlab as "datamat"
% nsubjects by nedges single

datamat_raw = load([datapath, 'train_newMSNnetwork.mat']).MSNnetwork; %% test
datamat = reshape(datamat_raw, [], size(datamat_raw, 3))';

mask=1:length(datamat);
mask_mat=reshape(mask,[sqrt(length(datamat)),sqrt(length(datamat))]);
grot=triu(mask_mat,1);
grot=grot(grot>0);
datamat=datamat(:,grot);

%% Prep behaviour file
% Add symp data as behav-file (CBCL)
%all_data = readtable([behaviorpath 'all_data.csv']);
train_data = readtable([behaviorpath 'train_data.csv']); %% test
all_data = readtable([behaviorpath 'all_data.csv']); %% test
train_subject_ids = train_data.Identifiers;
match_indices = ismember(all_data.Identifiers, train_subject_ids);
all_data_matched = all_data(match_indices, :);

% 方法1：分步条件替换（新手易理解）
% 初始化数值列（和 sex 列长度一致）
sex_numeric = zeros(height(all_data_matched), 1);  % 先默认设为 0（对应 Female）
% 替换 Male 为 1
sex_numeric(strcmp(all_data_matched.sex, 'Male')) = 1;
% 替换 Other 为 3
sex_numeric(strcmp(all_data_matched.sex, 'Other')) = 3;

% 将替换后的数值赋值回原 sex 列（或新建列，如 data.sex_num = sex_numeric）
all_data_matched.sex = sex_numeric;

% 将站点信息转为数字符号
site_numeric = zeros(height(all_data_matched), 1); % HBNsiteCBIC = 0
site_numeric(strcmp(all_data_matched.study_site, 'HBNsiteCUNY')) = 1;
site_numeric(strcmp(all_data_matched.study_site, 'HBNsiteRU')) = 2;
site_numeric(strcmp(all_data_matched.study_site, 'HBNsiteSI')) = 3;
all_data_matched.study_site = site_numeric;

CBCL = all_data_matched(: ,369:487);
symp_mat = table2array(CBCL);

covars = table2array(all_data_matched(:, [2:5, 368]));

symp_mat_scaled=zscore(symp_mat);

symp_names = CBCL.Properties.VariableNames;

writematrix(symp_mat, [plspath 'symp_mat_no_nans_raw.txt'])
writematrix(symp_mat_scaled, [plspath 'symp_mat_no_nans_scaled.txt'])
writecell(symp_names, [plspath 'symptom_names.txt'])

%% Deconfound age,sex, scanner, tsnr, fd_mean
addpath /Users/lizheng/Matlab_Plugin/FSLNets

covars(:,2)=categorical(covars(:,2)); %sex
covars(:,1)=zscore(covars(:,1)); %age
covars(:,3)=zscore(covars(:,3)); %elous
covars(:,4)=categorical(covars(:,4)); %site
covars(:,5)=zscore(covars(:,5)); % TIV


datamat_unconfound=nets_unconfound(datamat,covars);
datamat = single(datamat_unconfound);

save([plspath 'edges_deconfounded_matched_with_symp_STRUCTdata.mat'], 'datamat')

%% Create input session-file to bypass Nifti requirement of PLS toolbox
ids = all_data_matched.Identifiers';

bad_coords = [];
behavdata = [];
behavname = {};
coords = 1:size(datamat,2); %probably coords inside nifti file for selections of brain voxels? in that case, for me it will be 1:n_edges
create_datamat_info = struct('brain_mask_file', '', ...
    'normalize_volume_mean', 0);
create_ver = '2';
datafile = 'edges_deconfounded_matched_with_symp_STRUCTdata.mat';
dims = [100,100,100,100];
origin = [0,0,0];
selected_subjects = ones(1,size(datamat,1)); %list of ones for n_subjects

session_info = struct('description', 'hbn edges' , ...
     'pls_data_path', '', ...
      'dataset_path', '', ...
      'num_behavior', 0, ...
         'behavdata', behavdata, ...
         'behavname', {behavname}, ...
    'datamat_prefix', 'edges_INSERT_effects', ...
    'num_conditions', 1, ...
         'condition', {{'edges'}}, ...
      'num_subjects', size(datamat,1), ...
           'subject', {strcat(ids, '_*.*')}, ... %insert cell array of participant IDs_*.*
         'subj_name', {strcat(ids, '_')}, ... %insert cell array of participant IDs_
       'cond_filter', {{'*disco.nii'}}, ...
        'subj_files', {ids}, ... %insert cell array of paste0(participant IDs_cond_filter) i.e. nifti file
           'img_ext', '*.nii', ...
     'num_subj_init', -1);

singeprecision = 1;
voxel_size = [2,2,2]; 

save([plspath 'edges_deconfounded_matched_with_symp_STRUCTsessiondata.mat'], ...
    'behavdata','create_datamat_info','dims','session_info','behavname','create_ver','origin',...
    'singeprecision','bad_coords','coords','datafile','selected_subjects','voxel_size')

%% Open PLS gui to generate analysis file
addpath(genpath('/Users/lizheng/Matlab_Plugin/PLS-main'));
plsgui

%% Run PLS batch
% NB! Remember to manually change correlation mode to 8 in
% STRUCTanalysis.txt-file

batch_plsgui([plspath 'edges_INSERT_effects_STRUCTanalysis.txt'])

