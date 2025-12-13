%%% Read PLS output and save to matrix for plotting in R
%%% IV 2023
cd /Users/lizheng/Desktop/同步文件夹/博士研究课题/OHBM会议数据分析/Version2/output
% Import data
load('edges_INSERT_effects_STRUCTresult.mat');
%% Check significance of permutation test
permp = result.perm_result.sprob;
corr = result.lvcorrs; 
varexp = result.s.^2/sum(result.s.^2)*100;

siglv = table(corr, permp, varexp);
writetable(siglv, 'symptom_output/sig_sympLVs.csv')

vsc = zscore(result.vsc); % Symptom subject weights
usc = zscore(result.usc); % Connectivity subject weights
sc = table(vsc, usc);
writetable(sc, 'symptom_output/subject_weights.csv')
%% Loop through all sig LVs
for x = 1:length(permp)
    if permp(x) >= 0.05
        continue
    else
        edgeweights = result.u(:,x); % Connectivity weights each latent variable (LV)
        behavweights = result.v(:,x); % Symptom weights on each latent variable (LV)

        behcorrs=result.boot_result.orig_corr(:,x);
        behcorrs_ul=result.boot_result.ulcorr(:,x);
        behcorrs_ll=result.boot_result.llcorr(:,x);

        % Check bootstrap for significant edges
        % Threshold at pseudo-z/BSR > 3
        % Filter to obtain sig edges
        % Mask non-sig by making them 0

        pseudoz = result.boot_result.compare_u(:,x);

        % Flip LV1, LV3, LV4, and LV6 and LV7
        if x == 1 | x == 3 | x == 4 | x == 6
            edgeweights = edgeweights*-1;
            behavweights = behavweights*-1;
            behcorrs = behcorrs*-1;
            behcorrs_ul = behcorrs_ul*-1;
            behcorrs_ll = behcorrs_ll*-1;
            pseudoz = pseudoz*-1;
        end

        edgedata_raw = [edgeweights, pseudoz];

        pseudoz_sig = pseudoz;
        pseudoz_sig(abs(pseudoz) < 3) = 0;
        edgeweights_sig = edgeweights;
        edgeweights_sig(abs(pseudoz) < 3) = 0;

        edgedata_sig = [edgeweights_sig, pseudoz_sig];

        writematrix(edgedata_sig, ['symptom_output/sig_edges_sympLV', num2str(x), '_sigonly.csv'])
        writematrix(edgedata_raw, ['symptom_output/sig_edges_sympLV', num2str(x), '_raw.csv'])

        % Prep edge weights output
        [edge_ordered, i] = sort(edgeweights(:,1), 'descend', 'ComparisonMethod','abs');
        
        labels = readtable("she100.txt", ...
    'FileType', 'text', ...
    'ReadVariableNames', false, ...
    'TextType', 'string', ...
    'Delimiter', 'none'); 

        edgemat = triu(ones(100),1);
        edgemat_idx = triu(ones(100),1);

        edgemat_idx(edgemat > 0) = 1:4950;
        edgemat(edgemat > 0) = edgeweights_sig;

        edgedat = (edgemat + edgemat')./(eye(100)+1);

        for j = 1:length(i)
            [label_r(j), label_c(j)] = find(edgemat_idx == i(j));
        end

        labels_row = labels.Var1(label_r);
        labels_col = labels.Var1(label_c);

        edge_ordered_masked = edgeweights_sig(i,1);

        edge = table(edge_ordered_masked, labels_row, labels_col);

        writematrix(edgedat, ['symptom_output/edge_weights_matrixformat_symp_LV' num2str(x), '.csv'])
        writetable(edge, ['symptom_output/edge_weights_ordered_sig_sympLV', num2str(x), '.csv'])
        
        % Prep behaviour weights output
        [behav_ordered, i] = sort(behavweights(:,1), 'descend', 'ComparisonMethod','abs');

        % 1. 检测文件导入选项
        opts = detectImportOptions("sympnames_for_edges.txt");
        
        % 2. 关键设置：取消所有分隔符（核心解决拆分问题）
        opts.Delimiter = '';  % 空字符串表示无分隔符，整行作为一个单元格
        
        % 3. 确保只保留一列（避免生成空列）
        %opts.SelectedVariableNames = opts.VariableNames{1};  % 仅保留第一列

        sympnames = readcell("sympnames_for_edges.txt", opts);

        cbcl = readtable('cbcl_new.xlsx', 'ReadVariableNames', 1);

        %sympnames = extractAfter(sympnames(:,1),'CBCL_');
        sympnames_ordered = sympnames(i,1);

        [idx, locb]=ismember(sympnames_ordered(:,1),cbcl.Variable);

        behav = table(behav_ordered, sympnames_ordered);
        behav.itemtext = cbcl.Question(locb);

        writetable(behav, ['symptom_output/symptom_weights_ordered_sig_sympLV', num2str(x), '.csv'])

        % Prep behaviour loadings output
        [behcorrs_ordered, i] = sort(behcorrs(:,1), 'descend', 'ComparisonMethod','abs');

        behcorrs_ul_ordered = behcorrs_ul(i,1);
        behcorrs_ll_ordered = behcorrs_ll(i,1);

                % 1. 检测文件导入选项
        opts = detectImportOptions("sympnames_for_edges.txt");
        
        % 2. 关键设置：取消所有分隔符（核心解决拆分问题）
        opts.Delimiter = '';  % 空字符串表示无分隔符，整行作为一个单元格
        
        % 3. 确保只保留一列（避免生成空列）
        %opts.SelectedVariableNames = opts.VariableNames{1};  % 仅保留第一列

        sympnames = readcell("sympnames_for_edges.txt", opts);

        % sympnames = extractAfter(sympnames(:,1),'CBCL_');
        sympnames_ordered2 = sympnames(i,1);

        [idx, locb]=ismember(sympnames_ordered2(:,1),cbcl.Variable);

        behavcorrs = table(behcorrs_ordered, behcorrs_ul_ordered, behcorrs_ll_ordered, sympnames_ordered2);
        behavcorrs.itemtext = cbcl.Question(locb);

        writetable(behavcorrs, ['symptom_output/symptom_loadings_ordered_sig_sympLV', num2str(x), '.csv'])

    end
end