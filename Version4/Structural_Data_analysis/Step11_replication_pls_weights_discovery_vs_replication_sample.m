%%% Replication of PLS weights in discovery and replication sample
%%%
%%% Derive weights in replication sample based on PLS in discovery sample
%%% Correlate these with weights from PLS in replicaiton sample
%%% 
%%% Derive weights in discovery sample based on PLS in replication sample
%%% Correlate these with weights from PLS in discovery sample
%%%
%%% IV 2023
clear,clc
discpath='/Users/lizheng/Desktop/同步文件夹/博士研究课题/OHBM会议数据分析/Version4/output/';
replpath='/Users/lizheng/Desktop/同步文件夹/博士研究课题/OHBM会议数据分析/Version4/output/testdata/';
outpath='/Users/lizheng/Desktop/同步文件夹/博士研究课题/OHBM会议数据分析/Version4/validation/';

%% Load data
% Discovery sample
load([discpath 'edges_INSERT_effects_STRUCTresult.mat']);
disc_brainweights = result.u;
disc_behavweights = result.v;
disc_brain_subj_weights = result.usc;
disc_behav_subj_weights = result.vsc;

load([discpath 'edges_deconfounded_matched_with_symp_STRUCTdata.mat']);
disc_brain = datamat;
disc_behav = readmatrix([discpath 'symp_mat_no_nans_scaled.txt']);

% Replication sample
load([replpath 'edges_INSERT_effects_STRUCTresult.mat']);
repl_brainweights = result.u;
repl_behavweights = result.v;
repl_brain_subj_weights = result.usc;
repl_behav_subj_weights = result.vsc;

load([replpath 'edges_deconfounded_matched_with_symp_STRUCTdata.mat']);
repl_brain = datamat;
repl_behav = readmatrix([replpath 'symp_mat_no_nans_raw.txt']);

%% Loop correlation between brain weights

% Repl in disc
for i = 1:6
    disc = disc_brain_subj_weights(:,i);
    repl = disc_brain * repl_brainweights(:,i);
    
    [r(i), p(i)] = corr(disc, repl,'rows', 'pairwise'); 
    
    disc_out(:,i) = disc;
    repl_out(:,i) = repl; 

    % Permute
    grotALL=[];
    nperm=100;
    for j=1:nperm
        permcorr=corr(disc,repl(randperm(size(repl,1)),1), 'Type', 'Pearson');
        grotALL=[grotALL;permcorr'];
        disp(['Perm: ', num2str(j)])
    end

    nonperm_r = r(i);
    perm_r = grotALL;

    % Significance
    sig(i) = (sum(perm_r >= max(abs(nonperm_r)))+1) ./ nperm;
    
    histogram(perm_r); xline(max(abs(nonperm_r)),'--r',{'Non-permuted','R value'});
    saveas(gcf,[outpath filesep 'Repl_in_disc_brain_weights_corr_permresult_LV', num2str(i), '.png'])
    close

end

out = table(r', p', sig');
out.Properties.VariableNames = {'R', 'p', 'perm_p'};
writetable(out, [outpath filesep 'brain_weights_repl_in_disc_corr.csv'])

writematrix(disc_out,[outpath filesep 'brain_disc_subject_weights_repl_in_disc.csv'])
writematrix(repl_out,[outpath filesep 'brain_repl_subject_weights_repl_in_disc.csv'])

disc_out = [];
repl_out = [];

% Disc in repl
r = [];
p = [];

for i = 1:7
    disc = repl_brain * disc_brainweights(:,i);
    repl = repl_brain_subj_weights(:,i);
    
    [r(i), p(i)] = corr(disc, repl,'rows', 'pairwise'); 
    
    disc_out(:,i) = disc;
    repl_out(:,i) = repl; 
    
    % Permute
    grotALL=[];
    nperm=100;
    for j=1:nperm
        permcorr=corr(disc,repl(randperm(size(repl,1)),1), 'Type', 'Pearson');
        grotALL=[grotALL;permcorr'];
        disp(['Perm: ', num2str(j)])
    end

    nonperm_r = r(i);
    perm_r = grotALL;

    % Significance
    sig(i) = (sum(perm_r >= max(abs(nonperm_r)))+1) ./ nperm;

    histogram(perm_r); xline(max(abs(nonperm_r)),'--r',{'Non-permuted','R value'});
    saveas(gcf,[outpath filesep 'Disc_in_repl_brain_weights_corr_permresult_LV', num2str(i), '.png'])
    close

end

out = table(r', p', sig');
out.Properties.VariableNames = {'R', 'p', 'perm_p'};
writetable(out, [outpath filesep 'brain_weights_disc_in_repl_corr.csv'])

writematrix(disc_out,[outpath filesep 'brain_disc_subject_weights_disc_in_repl.csv'])
writematrix(repl_out,[outpath filesep 'brain_repl_subject_weights_disc_in_repl.csv'])

disc_out = [];
repl_out = [];

%% Loop correlation between behaviour/symptom weights

% Repl in disc
for i = 1:7
    disc = disc_behav_subj_weights(:,i);
    repl = disc_behav * repl_behavweights(:,i);
    
    [r(i), p(i)] = corr(disc, repl,'rows', 'pairwise'); 

    disc_out(:,i) = disc;
    repl_out(:,i) = repl; 

    % Permute
    grotALL=[];
    nperm=100;
    for j=1:nperm
        permcorr=corr(disc,repl(randperm(size(repl,1)),1), 'Type', 'Pearson');
        grotALL=[grotALL;permcorr'];
        disp(['Perm: ', num2str(j)])
    end

    nonperm_r = r(i);
    perm_r = grotALL;

    % Significance
    sig(i) = (sum(perm_r >= max(abs(nonperm_r)))+1) ./ nperm;
    
    histogram(perm_r); xline(max(abs(nonperm_r)),'--r',{'Non-permuted','R value'});
    saveas(gcf,[outpath filesep 'Repl_in_disc_behav_weights_corr_permresult_LV', num2str(i), '.png'])
    close

end

out = table(r', p', sig');
out.Properties.VariableNames = {'R', 'p', 'perm_p'};
writetable(out, [outpath filesep 'behav_weights_repl_in_disc_corr.csv'])

writematrix(disc_out,[outpath filesep 'behav_disc_subject_weights_repl_in_disc.csv'])
writematrix(repl_out,[outpath filesep 'behav_repl_subject_weights_repl_in_disc.csv'])

disc_out = [];
repl_out = [];

% Disc in repl
r = [];
p = [];

for i = 1:7
    disc = repl_behav * disc_behavweights(:,i);
    repl = repl_behav_subj_weights(:,i);
    
    [r(i), p(i)] = corr(disc, repl,'rows', 'pairwise'); 
    
    disc_out(:,i) = disc;
    repl_out(:,i) = repl; 

    % Permute
    grotALL=[];
    nperm=100;
    for j=1:nperm
        permcorr=corr(disc,repl(randperm(size(repl,1)),1), 'Type', 'Pearson');
        grotALL=[grotALL;permcorr'];
        disp(['Perm: ', num2str(j)])
    end

    nonperm_r = r(i);
    perm_r = grotALL;

    % Significance
    sig(i) = (sum(perm_r >= max(abs(nonperm_r)))+1) ./ nperm;

    histogram(perm_r); xline(max(abs(nonperm_r)),'--r',{'Non-permuted','R value'});
    saveas(gcf,[outpath filesep 'Disc_in_repl_behav_weights_corr_permresult_LV', num2str(i), '.png'])
    close

end

out = table(r', p', sig');
out.Properties.VariableNames = {'R', 'p', 'perm_p'};
writetable(out, [outpath filesep 'behav_weights_disc_in_repl_corr.csv'])

writematrix(disc_out,[outpath filesep 'behav_disc_subject_weights_disc_in_repl.csv'])
writematrix(repl_out,[outpath filesep 'behav_repl_subject_weights_disc_in_repl.csv'])
