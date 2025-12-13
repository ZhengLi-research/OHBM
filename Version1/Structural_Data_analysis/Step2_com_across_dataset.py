# combine all 5 datasets
# 对不同的脑结构指标，整合来自不同数据集的数据
# 还没有完全自动化，后续需要迭代
import pandas as pd

datapath = '/Users/lizheng/Desktop/code/RBCcode/rbc_data_analysis/Newdata/'


dtype = ['ct_pass', 'gv_pass', 'sa_pass', 'mc_pass', 'gc_pass', 'fi_pass', 'ci_pass', 'sv_pass']

for iType in dtype:
    pnc = pd.read_csv(datapath + '/pnc_df_%s.tsv' % iType,
                      delimiter='\t')
    hbn = pd.read_csv(datapath + '/hbn_df_%s.tsv' % iType,
                      delimiter='\t')
    # bhrc = pd.read_csv(datapath + 'data/dataR/bhrc_df_%s.tsv' % iType,
    #                    delimiter='\t')
    # nki = pd.read_csv(datapath + 'data/dataR/nki_df_%s.tsv' % iType,
    #                   delimiter='\t')
    # ccnp = pd.read_csv(datapath + 'data/dataR/ccnp_df_%s.tsv' % iType,
    #                    delimiter='\t')

    combined_df = pd.concat([pnc, hbn], axis=0) #bhrc, nki, ccnp

    combined_df.reset_index(inplace=True, drop=True)

    combined_df.to_csv(datapath + '/combined_pnc_hbn_df_%s.tsv' % iType,
                       index=False, sep='\t')