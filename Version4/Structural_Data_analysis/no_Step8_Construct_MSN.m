%Importing the data:
clear, clc,close
cd /Users/lizheng/Desktop/同步文件夹/博士研究课题/OHBM会议数据分析/Version2/Structural_data/Structural_split/

a=load('ci_test.mat').fin;
b=load('fi_test.mat').fin;
c=load('gc_test.mat').fin;
d=load('gv_test.mat').fin;
e=load('mc_test.mat').fin;
f=load('sa_test.mat').fin;
g=load('ct_test.mat').fin;

a=table2array(a);
b=table2array(b);
c=table2array(c);
d=table2array(d);
e=table2array(e);
f=table2array(f);
g=table2array(g);


nregs=210; % number of regions
nsubs=size(a,1); % calculate the number of subjects based 'CurvInd.mat'


% calculate Z-score for each freesurfer index z-score the inputs

PARC500_CT_zscore=zscore(transpose(a));
PARC500_SA_zscore=zscore(transpose(b));
PARC500_GM_zscore=zscore(transpose(c));
PARC500_MC_zscore=zscore(transpose(d));
PARC500_GC_zscore=zscore(transpose(e));
PARC500_FA_zscore=zscore(transpose(f));
PARC500_MD_zscore=zscore(transpose(g));

% Create a cell for each subject with all of the required inputs:
clear subj_features7
for subj=1:nsubs
    subj_features7{1,subj}(:,1)=PARC500_CT_zscore(:,subj);
    subj_features7{1,subj}(:,2)=PARC500_SA_zscore(:,subj);
    subj_features7{1,subj}(:,3)=PARC500_GM_zscore(:,subj);
    subj_features7{1,subj}(:,4)=PARC500_MC_zscore(:,subj);
    subj_features7{1,subj}(:,5)=PARC500_GC_zscore(:,subj);
    subj_features7{1,subj}(:,6)=PARC500_FA_zscore(:,subj);
    subj_features7{1,subj}(:,7)=PARC500_MD_zscore(:,subj);
end

% Calculate the MS matrices by correlating all inputs and set the diagonal to zero:
for subj=1:nsubs
    subj_MSN_7{1,subj}=corr(transpose(subj_features7{1,subj}));
    subj_MSN_7{1,subj}(logical(eye(size(subj_MSN_7{1,subj})))) = 0;
    Z_subj_MSN_7{1,subj}=0.5*log((1+subj_MSN_7{1,subj})./(1-subj_MSN_7{1,subj}));   
end

% clear meanMS_regional
% for subj=1:nsubs
%     meanMS_regional(subj,:)=sum(subj_MSN_7{1,subj})./(nregs-1);
% end

MSNnetwork=cat(3,Z_subj_MSN_7{:});
save('test_newMSNnetwork.mat',"MSNnetwork",'-mat')
