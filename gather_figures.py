from pathlib import Path

FIG2A = ['figures/bar/expression_os.tsv', 'figures/bar/expression_pfi.tsv']

FIG2B = ['figures/bar/microbial_os.tsv', 'figures/bar/microbial_pfi.tsv']

FIG2C = [
    'figures/violin/tcga_acc_surv_os_htseq_counts_cnet_violin.tsv',
    'figures/violin/tcga_cesc_surv_os_htseq_counts_cnet_violin.tsv',
    'figures/violin/tcga_esca_surv_os_htseq_counts_cnet_violin.tsv',
    'figures/violin/tcga_kirc_surv_os_htseq_counts_cnet_violin.tsv',
    'figures/violin/tcga_acc_surv_pfi_htseq_counts_cnet_violin.tsv',
    'figures/violin/tcga_lgg_surv_pfi_htseq_counts_cnet_violin.tsv',
    'figures/violin/tcga_acc_surv_os_kraken_cnet_violin.tsv',
    'figures/violin/tcga_cesc_surv_os_kraken_cnet_violin.tsv',
    'figures/violin/tcga_esca_surv_os_kraken_cnet_violin.tsv',
    'figures/violin/tcga_kirc_surv_os_kraken_cnet_violin.tsv',
    'figures/violin/tcga_acc_surv_pfi_kraken_cnet_violin.tsv',
    'figures/violin/tcga_lgg_surv_pfi_kraken_cnet_violin.tsv',
]

FIG3A = ['figures/bar/microbial_response.tsv']

FIG3B = ['figures/bar/kraken_balanced_accuracy_bar_comp.tsv']

FIG3C = [
    'figures/violin/tcga_brca_resp_docetaxel_kraken_limma_violin.tsv',
    'figures/violin/tcga_brca_resp_docetaxel_kraken_rfe_violin.tsv',
    'figures/violin/tcga_stad_resp_cisplatin_kraken_limma_violin.tsv',
    'figures/violin/tcga_stad_resp_leucovorin_kraken_limma_violin.tsv',
    'figures/violin/tcga_stad_resp_oxaliplatin_kraken_limma_violin.tsv',
]

FIG3D = [
    'figures/roc_pr/tcga_brca_resp_docetaxel_kraken_limma_roc_auc.tsv',
    'figures/roc_pr/tcga_brca_resp_docetaxel_kraken_rfe_roc_auc.tsv',
    'figures/roc_pr/tcga_stad_resp_cisplatin_kraken_limma_roc_auc.tsv',
    'figures/roc_pr/tcga_stad_resp_leucovorin_kraken_limma_roc_auc.tsv',
    'figures/roc_pr/tcga_stad_resp_oxaliplatin_kraken_limma_roc_auc.tsv',
]

FIG3E = [
    'figures/roc_pr/tcga_brca_resp_docetaxel_kraken_limma_pr_auc.tsv',
    'figures/roc_pr/tcga_brca_resp_docetaxel_kraken_rfe_pr_auc.tsv',
    'figures/roc_pr/tcga_stad_resp_cisplatin_kraken_limma_pr_auc.tsv',
    'figures/roc_pr/tcga_stad_resp_leucovorin_kraken_limma_pr_auc.tsv',
    'figures/roc_pr/tcga_stad_resp_oxaliplatin_kraken_limma_pr_auc.tsv',
]

FIG4A = ['figures/bar/expression_response.tsv']

FIG4B = ['figures/bar/htseq_balanced_accuracy_bar_comp.tsv']

FIG4C = [
    'figures/violin/tcga_blca_resp_cisplatin_htseq_counts_rfe_violin.tsv',
    'figures/violin/tcga_blca_resp_gemcitabine_htseq_counts_rfe_violin.tsv',
    'figures/violin/tcga_hnsc_resp_carboplatin_htseq_counts_rfe_violin.tsv',
    'figures/violin/tcga_sarc_resp_docetaxel_htseq_counts_edger_violin.tsv',
    'figures/violin/tcga_sarc_resp_gemcitabine_htseq_counts_rfe_violin.tsv',
    'figures/violin/tcga_tgct_resp_bleomycin_htseq_counts_lgr_violin.tsv',
]

FIG4D = [
    'figures/roc_pr/tcga_blca_resp_cisplatin_htseq_counts_rfe_roc_auc.tsv',
    'figures/roc_pr/tcga_blca_resp_gemcitabine_htseq_counts_rfe_roc_auc.tsv',
    'figures/roc_pr/tcga_hnsc_resp_carboplatin_htseq_counts_rfe_roc_auc.tsv',
    'figures/roc_pr/tcga_sarc_resp_docetaxel_htseq_counts_edger_roc_auc.tsv',
    'figures/roc_pr/tcga_sarc_resp_gemcitabine_htseq_counts_rfe_roc_auc.tsv',
    'figures/roc_pr/tcga_tgct_resp_bleomycin_htseq_counts_lgr_roc_auc.tsv',
]

FIG4E = [
    'figures/roc_pr/tcga_blca_resp_cisplatin_htseq_counts_rfe_pr_auc.tsv',
    'figures/roc_pr/tcga_blca_resp_gemcitabine_htseq_counts_rfe_pr_auc.tsv',
    'figures/roc_pr/tcga_hnsc_resp_carboplatin_htseq_counts_rfe_pr_auc.tsv',
    'figures/roc_pr/tcga_sarc_resp_docetaxel_htseq_counts_edger_pr_auc.tsv',
    'figures/roc_pr/tcga_sarc_resp_gemcitabine_htseq_counts_rfe_pr_auc.tsv',
    'figures/roc_pr/tcga_tgct_resp_bleomycin_htseq_counts_lgr_pr_auc.tsv',
]

FIG5AC = 'figures/mlcomp/venn_summary.tsv'

FIG5BD = 'figures/mlcomp/feature_correlation.tsv'

FIG6A = [
    'figures/perm_hist/tcga_brca_resp_docetaxel_kraken_limma_perm_hist.tsv',
    'figures/perm_hist/tcga_sarc_resp_docetaxel_kraken_rfe_perm_hist.tsv',
    'figures/perm_hist/tcga_stad_resp_cisplatin_kraken_limma_perm_hist.tsv',
    'figures/perm_hist/tcga_stad_resp_leucovorin_kraken_limma_perm_hist.tsv',
    'figures/perm_hist/tcga_stad_resp_oxaliplatin_kraken_limma_perm_hist.tsv',
]

FIG6B = [
    'figures/roc_pr/tcga_brca_resp_docetaxel_kraken_limma_k_vs_score.tsv',
    'figures/roc_pr/tcga_sarc_resp_docetaxel_kraken_rfe_k_vs_score.tsv',
    'figures/roc_pr/tcga_stad_resp_cisplatin_kraken_limma_k_vs_score.tsv',
    'figures/roc_pr/tcga_stad_resp_leucovorin_kraken_limma_k_vs_score.tsv',
    'figures/roc_pr/tcga_stad_resp_oxaliplatin_kraken_limma_k_vs_score.tsv',
]

FIG6C = [
    'figures/perm_hist/tcga_blca_resp_cisplatin_htseq_counts_rfe_perm_hist.tsv',
    'figures/perm_hist/tcga_blca_resp_gemcitabine_htseq_counts_rfe_perm_hist.tsv',
    'figures/perm_hist/tcga_hnsc_resp_carboplatin_htseq_counts_rfe_perm_hist.tsv',
    'figures/perm_hist/tcga_sarc_resp_docetaxel_htseq_counts_edger_perm_hist.tsv',
    'figures/perm_hist/tcga_sarc_resp_gemcitabine_htseq_counts_rfe_perm_hist.tsv',
    'figures/perm_hist/tcga_tgct_resp_bleomycin_htseq_counts_lgr_perm_hist.tsv',
]

FIG6D = [
    'figures/roc_pr/tcga_blca_resp_cisplatin_htseq_counts_rfe_k_vs_score.tsv',
    'figures/roc_pr/tcga_blca_resp_gemcitabine_htseq_counts_rfe_k_vs_score.tsv',
    'figures/roc_pr/tcga_hnsc_resp_carboplatin_htseq_counts_rfe_k_vs_score.tsv',
    'figures/roc_pr/tcga_sarc_resp_docetaxel_htseq_counts_edger_k_vs_score.tsv',
    'figures/roc_pr/tcga_sarc_resp_gemcitabine_htseq_counts_rfe_k_vs_score.tsv',
    'figures/roc_pr/tcga_tgct_resp_bleomycin_htseq_counts_lgr_c_vs_score.tsv',
]

FIGS1A = [
    'figures/violin/tcga_acc_surv_os_htseq_counts_cnet_violin.tsv',
    'figures/violin/tcga_blca_surv_os_htseq_counts_cnet_violin.tsv',
    'figures/violin/tcga_cesc_surv_os_htseq_counts_cnet_violin.tsv',
    'figures/violin/tcga_hnsc_surv_os_htseq_counts_cnet_violin.tsv',
    'figures/violin/tcga_kirc_surv_os_htseq_counts_cnet_violin.tsv',
    'figures/violin/tcga_kirp_surv_os_htseq_counts_cnet_violin.tsv',
    'figures/violin/tcga_lgg_surv_os_htseq_counts_cnet_violin.tsv',
    'figures/violin/tcga_lihc_surv_os_htseq_counts_cnet_violin.tsv',
    'figures/violin/tcga_luad_surv_os_htseq_counts_cnet_violin.tsv',
    'figures/violin/tcga_meso_surv_os_htseq_counts_cnet_violin.tsv',
    'figures/violin/tcga_pcpg_surv_os_htseq_counts_cnet_violin.tsv',
    'figures/violin/tcga_sarc_surv_os_htseq_counts_cnet_violin.tsv',
    'figures/violin/tcga_skcm_surv_os_htseq_counts_cnet_violin.tsv',
    'figures/violin/tcga_thym_surv_os_htseq_counts_cnet_violin.tsv',
    'figures/violin/tcga_ucec_surv_os_htseq_counts_cnet_violin.tsv',
    'figures/violin/tcga_uvm_surv_os_htseq_counts_cnet_violin.tsv',
]

FIGS1B = [
    'figures/td_auc/tcga_acc_surv_os_htseq_counts_cnet_td_auc.tsv',
    'figures/td_auc/tcga_blca_surv_os_htseq_counts_cnet_td_auc.tsv',
    'figures/td_auc/tcga_cesc_surv_os_htseq_counts_cnet_td_auc.tsv',
    'figures/td_auc/tcga_hnsc_surv_os_htseq_counts_cnet_td_auc.tsv',
    'figures/td_auc/tcga_kirc_surv_os_htseq_counts_cnet_td_auc.tsv',
    'figures/td_auc/tcga_kirp_surv_os_htseq_counts_cnet_td_auc.tsv',
    'figures/td_auc/tcga_lgg_surv_os_htseq_counts_cnet_td_auc.tsv',
    'figures/td_auc/tcga_lihc_surv_os_htseq_counts_cnet_td_auc.tsv',
    'figures/td_auc/tcga_luad_surv_os_htseq_counts_cnet_td_auc.tsv',
    'figures/td_auc/tcga_meso_surv_os_htseq_counts_cnet_td_auc.tsv',
    'figures/td_auc/tcga_pcpg_surv_os_htseq_counts_cnet_td_auc.tsv',
    'figures/td_auc/tcga_sarc_surv_os_htseq_counts_cnet_td_auc.tsv',
    'figures/td_auc/tcga_skcm_surv_os_htseq_counts_cnet_td_auc.tsv',
    'figures/td_auc/tcga_thym_surv_os_htseq_counts_cnet_td_auc.tsv',
    'figures/td_auc/tcga_ucec_surv_os_htseq_counts_cnet_td_auc.tsv',
    'figures/td_auc/tcga_uvm_surv_os_htseq_counts_cnet_td_auc.tsv',
]

FIGS2A = [
    'figures/violin/tcga_acc_surv_pfi_htseq_counts_cnet_violin.tsv',
    'figures/violin/tcga_blca_surv_pfi_htseq_counts_cnet_violin.tsv',
    'figures/violin/tcga_cesc_surv_pfi_htseq_counts_cnet_violin.tsv',
    'figures/violin/tcga_kirc_surv_pfi_htseq_counts_cnet_violin.tsv',
    'figures/violin/tcga_kirp_surv_pfi_htseq_counts_cnet_violin.tsv',
    'figures/violin/tcga_lgg_surv_pfi_htseq_counts_cnet_violin.tsv',
    'figures/violin/tcga_lihc_surv_pfi_htseq_counts_cnet_violin.tsv',
    'figures/violin/tcga_lusc_surv_pfi_htseq_counts_cnet_violin.tsv',
    'figures/violin/tcga_meso_surv_pfi_htseq_counts_cnet_violin.tsv',
    'figures/violin/tcga_paad_surv_pfi_htseq_counts_cnet_violin.tsv',
    'figures/violin/tcga_prad_surv_pfi_htseq_counts_cnet_violin.tsv',
    'figures/violin/tcga_sarc_surv_pfi_htseq_counts_cnet_violin.tsv',
    'figures/violin/tcga_skcm_surv_pfi_htseq_counts_cnet_violin.tsv',
    'figures/violin/tcga_stad_surv_pfi_htseq_counts_cnet_violin.tsv',
    'figures/violin/tcga_ucec_surv_pfi_htseq_counts_cnet_violin.tsv',
    'figures/violin/tcga_uvm_surv_pfi_htseq_counts_cnet_violin.tsv',
]

FIGS2B = [
    'figures/td_auc/tcga_acc_surv_pfi_htseq_counts_cnet_td_auc.tsv',
    'figures/td_auc/tcga_blca_surv_pfi_htseq_counts_cnet_td_auc.tsv',
    'figures/td_auc/tcga_cesc_surv_pfi_htseq_counts_cnet_td_auc.tsv',
    'figures/td_auc/tcga_kirc_surv_pfi_htseq_counts_cnet_td_auc.tsv',
    'figures/td_auc/tcga_kirp_surv_pfi_htseq_counts_cnet_td_auc.tsv',
    'figures/td_auc/tcga_lgg_surv_pfi_htseq_counts_cnet_td_auc.tsv',
    'figures/td_auc/tcga_lihc_surv_pfi_htseq_counts_cnet_td_auc.tsv',
    'figures/td_auc/tcga_lusc_surv_pfi_htseq_counts_cnet_td_auc.tsv',
    'figures/td_auc/tcga_meso_surv_pfi_htseq_counts_cnet_td_auc.tsv',
    'figures/td_auc/tcga_paad_surv_pfi_htseq_counts_cnet_td_auc.tsv',
    'figures/td_auc/tcga_prad_surv_pfi_htseq_counts_cnet_td_auc.tsv',
    'figures/td_auc/tcga_sarc_surv_pfi_htseq_counts_cnet_td_auc.tsv',
    'figures/td_auc/tcga_skcm_surv_pfi_htseq_counts_cnet_td_auc.tsv',
    'figures/td_auc/tcga_stad_surv_pfi_htseq_counts_cnet_td_auc.tsv',
    'figures/td_auc/tcga_ucec_surv_pfi_htseq_counts_cnet_td_auc.tsv',
    'figures/td_auc/tcga_uvm_surv_pfi_htseq_counts_cnet_td_auc.tsv',
]

FIGS3A = [
    'figures/violin/tcga_acc_surv_os_kraken_cnet_violin.tsv',
    'figures/violin/tcga_cesc_surv_os_kraken_cnet_violin.tsv',
    'figures/violin/tcga_esca_surv_os_kraken_cnet_violin.tsv',
    'figures/violin/tcga_kirc_surv_os_kraken_cnet_violin.tsv',
    'figures/violin/tcga_acc_surv_pfi_kraken_cnet_violin.tsv',
    'figures/violin/tcga_lgg_surv_pfi_kraken_cnet_violin.tsv',
]

FIGS3B = [
    'figures/td_auc/tcga_acc_surv_os_kraken_cnet_td_auc.tsv',
    'figures/td_auc/tcga_cesc_surv_os_kraken_cnet_td_auc.tsv',
    'figures/td_auc/tcga_esca_surv_os_kraken_cnet_td_auc.tsv',
    'figures/td_auc/tcga_kirc_surv_os_kraken_cnet_td_auc.tsv',
    'figures/td_auc/tcga_acc_surv_pfi_kraken_cnet_td_auc.tsv',
    'figures/td_auc/tcga_lgg_surv_pfi_kraken_cnet_td_auc.tsv',
]

FIGS4AI = [
    'figures/violin/tcga_sarc_surv_os_combo_cnet_violin_comp.tsv',
    'figures/violin/tcga_stad_surv_pfi_combo_cnet_violin_comp.tsv',
    'figures/violin/tcga_thym_surv_os_combo_cnet_violin_comp.tsv',
]

FIGS4AII = [
    'figures/violin/tcga_sarc_surv_os_combo_cnet_violin.tsv',
    'figures/violin/tcga_stad_surv_pfi_combo_cnet_violin.tsv',
    'figures/violin/tcga_thym_surv_os_combo_cnet_violin.tsv',
]

FIGS4B = ['figures/bar/htseq_balanced_accuracy_bar_comp.tsv']

FIGS4CI = ['figures/violin/tcga_blca_resp_cisplatin_combo_rfe_violin_comp.tsv']
FIGS4CII = ['figures/violin/tcga_blca_resp_cisplatin_combo_rfe_violin.tsv']

FIGS5A = [
    'figures/perm_hist/tcga_brca_resp_docetaxel_kraken_lgr_perm_hist.tsv',
    'figures/perm_hist/tcga_brca_resp_docetaxel_kraken_limma_perm_hist.tsv',
    'figures/perm_hist/tcga_brca_resp_docetaxel_kraken_rfe_perm_hist.tsv',
    'figures/perm_hist/tcga_sarc_resp_docetaxel_kraken_lgr_perm_hist.tsv',
    'figures/perm_hist/tcga_sarc_resp_docetaxel_kraken_limma_perm_hist.tsv',
    'figures/perm_hist/tcga_sarc_resp_docetaxel_kraken_rfe_perm_hist.tsv',
    'figures/perm_hist/tcga_stad_resp_cisplatin_kraken_lgr_perm_hist.tsv',
    'figures/perm_hist/tcga_stad_resp_cisplatin_kraken_limma_perm_hist.tsv',
    'figures/perm_hist/tcga_stad_resp_cisplatin_kraken_rfe_perm_hist.tsv',
    'figures/perm_hist/tcga_stad_resp_leucovorin_kraken_lgr_perm_hist.tsv',
    'figures/perm_hist/tcga_stad_resp_leucovorin_kraken_limma_perm_hist.tsv',
    'figures/perm_hist/tcga_stad_resp_leucovorin_kraken_rfe_perm_hist.tsv',
    'figures/perm_hist/tcga_stad_resp_oxaliplatin_kraken_lgr_perm_hist.tsv',
    'figures/perm_hist/tcga_stad_resp_oxaliplatin_kraken_limma_perm_hist.tsv',
    'figures/perm_hist/tcga_stad_resp_oxaliplatin_kraken_rfe_perm_hist.tsv',
]

FIGS5B = [
    'figures/roc_pr/tcga_brca_resp_docetaxel_kraken_lgr_c_vs_score.tsv',
    'figures/roc_pr/tcga_brca_resp_docetaxel_kraken_lgr_l1r_vs_score.tsv',
    'figures/roc_pr/tcga_brca_resp_docetaxel_kraken_limma_k_vs_score.tsv',
    'figures/roc_pr/tcga_brca_resp_docetaxel_kraken_rfe_k_vs_score.tsv',
    'figures/roc_pr/tcga_sarc_resp_docetaxel_kraken_lgr_c_vs_score.tsv',
    'figures/roc_pr/tcga_sarc_resp_docetaxel_kraken_lgr_l1r_vs_score.tsv',
    'figures/roc_pr/tcga_sarc_resp_docetaxel_kraken_limma_k_vs_score.tsv',
    'figures/roc_pr/tcga_sarc_resp_docetaxel_kraken_rfe_k_vs_score.tsv',
    'figures/roc_pr/tcga_stad_resp_cisplatin_kraken_lgr_c_vs_score.tsv',
    'figures/roc_pr/tcga_stad_resp_cisplatin_kraken_lgr_l1r_vs_score.tsv',
    'figures/roc_pr/tcga_stad_resp_cisplatin_kraken_limma_k_vs_score.tsv',
    'figures/roc_pr/tcga_stad_resp_cisplatin_kraken_rfe_k_vs_score.tsv',
    'figures/roc_pr/tcga_stad_resp_leucovorin_kraken_lgr_c_vs_score.tsv',
    'figures/roc_pr/tcga_stad_resp_leucovorin_kraken_lgr_l1r_vs_score.tsv',
    'figures/roc_pr/tcga_stad_resp_leucovorin_kraken_limma_k_vs_score.tsv',
    'figures/roc_pr/tcga_stad_resp_leucovorin_kraken_rfe_k_vs_score.tsv',
    'figures/roc_pr/tcga_stad_resp_oxaliplatin_kraken_lgr_c_vs_score.tsv',
    'figures/roc_pr/tcga_stad_resp_oxaliplatin_kraken_lgr_l1r_vs_score.tsv',
    'figures/roc_pr/tcga_stad_resp_oxaliplatin_kraken_limma_k_vs_score.tsv',
    'figures/roc_pr/tcga_stad_resp_oxaliplatin_kraken_rfe_k_vs_score.tsv',
]

FIGS6A = [
    'figures/perm_hist/tcga_blca_resp_cisplatin_htseq_counts_rfe_perm_hist.tsv',
    'figures/perm_hist/tcga_blca_resp_cisplatin_htseq_counts_edger_perm_hist.tsv',
    'figures/perm_hist/tcga_blca_resp_cisplatin_htseq_counts_lgr_perm_hist.tsv',
    'figures/perm_hist/tcga_blca_resp_gemcitabine_htseq_counts_rfe_perm_hist.tsv',
    'figures/perm_hist/tcga_blca_resp_gemcitabine_htseq_counts_edger_perm_hist.tsv',
    'figures/perm_hist/tcga_blca_resp_gemcitabine_htseq_counts_lgr_perm_hist.tsv',
    'figures/perm_hist/tcga_hnsc_resp_carboplatin_htseq_counts_rfe_perm_hist.tsv',
    'figures/perm_hist/tcga_hnsc_resp_carboplatin_htseq_counts_edger_perm_hist.tsv',
    'figures/perm_hist/tcga_hnsc_resp_carboplatin_htseq_counts_lgr_perm_hist.tsv',
    'figures/perm_hist/tcga_sarc_resp_docetaxel_htseq_counts_rfe_perm_hist.tsv',
    'figures/perm_hist/tcga_sarc_resp_docetaxel_htseq_counts_edger_perm_hist.tsv',
    'figures/perm_hist/tcga_sarc_resp_docetaxel_htseq_counts_lgr_perm_hist.tsv',
    'figures/perm_hist/tcga_sarc_resp_gemcitabine_htseq_counts_rfe_perm_hist.tsv',
    'figures/perm_hist/tcga_sarc_resp_gemcitabine_htseq_counts_edger_perm_hist.tsv',
    'figures/perm_hist/tcga_sarc_resp_gemcitabine_htseq_counts_lgr_perm_hist.tsv',
    'figures/perm_hist/tcga_tgct_resp_bleomycin_htseq_counts_rfe_perm_hist.tsv',
    'figures/perm_hist/tcga_tgct_resp_bleomycin_htseq_counts_edger_perm_hist.tsv',
    'figures/perm_hist/tcga_tgct_resp_bleomycin_htseq_counts_lgr_perm_hist.tsv',
]

FIGS6B = [
    'figures/roc_pr/tcga_blca_resp_cisplatin_htseq_counts_rfe_k_vs_score.tsv',
    'figures/roc_pr/tcga_blca_resp_cisplatin_htseq_counts_edger_k_vs_score.tsv',
    'figures/roc_pr/tcga_blca_resp_cisplatin_htseq_counts_lgr_c_vs_score.tsv',
    'figures/roc_pr/tcga_blca_resp_cisplatin_htseq_counts_lgr_l1r_vs_score.tsv',
    'figures/roc_pr/tcga_blca_resp_gemcitabine_htseq_counts_rfe_k_vs_score.tsv',
    'figures/roc_pr/tcga_blca_resp_gemcitabine_htseq_counts_edger_k_vs_score.tsv',
    'figures/roc_pr/tcga_blca_resp_gemcitabine_htseq_counts_lgr_c_vs_score.tsv',
    'figures/roc_pr/tcga_blca_resp_gemcitabine_htseq_counts_lgr_l1r_vs_score.tsv',
    'figures/roc_pr/tcga_hnsc_resp_carboplatin_htseq_counts_rfe_k_vs_score.tsv',
    'figures/roc_pr/tcga_hnsc_resp_carboplatin_htseq_counts_edger_k_vs_score.tsv',
    'figures/roc_pr/tcga_hnsc_resp_carboplatin_htseq_counts_lgr_c_vs_score.tsv',
    'figures/roc_pr/tcga_hnsc_resp_carboplatin_htseq_counts_lgr_l1r_vs_score.tsv',
    'figures/roc_pr/tcga_sarc_resp_docetaxel_htseq_counts_rfe_k_vs_score.tsv',
    'figures/roc_pr/tcga_sarc_resp_docetaxel_htseq_counts_edger_k_vs_score.tsv',
    'figures/roc_pr/tcga_sarc_resp_docetaxel_htseq_counts_lgr_c_vs_score.tsv',
    'figures/roc_pr/tcga_sarc_resp_docetaxel_htseq_counts_lgr_l1r_vs_score.tsv',
    'figures/roc_pr/tcga_sarc_resp_gemcitabine_htseq_counts_rfe_k_vs_score.tsv',
    'figures/roc_pr/tcga_sarc_resp_gemcitabine_htseq_counts_edger_k_vs_score.tsv',
    'figures/roc_pr/tcga_sarc_resp_gemcitabine_htseq_counts_lgr_c_vs_score.tsv',
    'figures/roc_pr/tcga_sarc_resp_gemcitabine_htseq_counts_lgr_l1r_vs_score.tsv',
    'figures/roc_pr/tcga_tgct_resp_bleomycin_htseq_counts_rfe_k_vs_score.tsv',
    'figures/roc_pr/tcga_tgct_resp_bleomycin_htseq_counts_edger_k_vs_score.tsv',
    'figures/roc_pr/tcga_tgct_resp_bleomycin_htseq_counts_lgr_c_vs_score.tsv',
    'figures/roc_pr/tcga_tgct_resp_bleomycin_htseq_counts_lgr_l1r_vs_score.tsv',
]

DATA_TYPE_NAMES = {
    'kraken': 'microbiome',
    'htseq': 'expression',
    'combo':   'combo',
}


def gather_separate_stats_fh(ofh, filenames):
    data_type_names = DATA_TYPE_NAMES
    for fileno, filename in enumerate(filenames):
        with open(filename) as lines:
            for iline, line in enumerate(lines):
                if fileno > 0 and iline == 0:
                    continue
                line = line.rstrip()
                F = line.split('\t')
                F[2] = data_type_names.get(F[2], F[2])
                print(*F, sep='\t', file=ofh)
    print("\n", file=ofh)

    stats_files = [f.replace('.tsv', '_stats.tsv') for f in filenames]
    for fileno, filename in enumerate(stats_files):
        with open(filename) as fh:
            for iline, line in enumerate(fh):
                if fileno > 0 and iline == 0:
                    continue
                line = line.rstrip()
                F = line.split('\t')
                F[0] = F[0].replace(' ', '_')
                print(*F, sep="\t", file=ofh)


def gather_separate_stats(outfile, filenames):
    with open(outfile, 'w') as ofh:
        gather_separate_stats_fh(ofh, filenames)


def gather_parse_filename_fh(ofh, filenames):
    data_type_names = DATA_TYPE_NAMES

    for fileno, filename in enumerate(filenames):
        with open(filename) as fh:
            pfilename = Path(filename)
            (_, cancer, undef, target, raw_data_type,
                    *model) = pfilename.stem.split('_')
            cancer = cancer.upper()
            model = model[0] if model[0] != 'counts' else model[1]
            data_type = data_type_names[raw_data_type]
            
            for iline, line in enumerate(fh):
                line = line.rstrip()
                F = line.split('\t')
                if iline == 0:
                    if fileno == 0:
                        print(
                            'cancer', 'target', 'data_type',
                            'model_code', *F, sep="\t", file=ofh)
                    continue
                print(cancer, target, data_type, model, *F,
                      sep="\t", file=ofh)

def gather_parse_filename(outfile, filenames):
    with open(outfile, 'w') as ofh:
        gather_parse_filename_fh(ofh, filenames)


def gather_parse_column_fh(ofh, filenames):
    data_type_names = DATA_TYPE_NAMES

    for fileno, filename in enumerate(filenames):
        with open(filename) as fh:
            for iline, line in enumerate(fh):
                line = line.rstrip()
                F = line.split('\t')
                if iline == 0:
                    if fileno == 0:
                        F = F[1:]
                        print(
                            'cancer', 'target', 'data_type',
                            'model_code', *F, sep="\t", file=ofh)
                    continue
                what_data = F[0]
                F = F[1:]
                (_, cancer, undef, target, raw_data_type,
                    *model) = what_data.split('_')
                cancer = cancer.upper()
                model = model[0] if model[0] != 'counts' else model[1]
                data_type = data_type_names[raw_data_type]
                print(cancer, target, data_type, model, *F,
                      sep="\t", file=ofh)


def gather_parse_column(outfile, filenames):
    with open(outfile, 'w') as ofh:
        gather_parse_column_fh(ofh, filenames)


def gather_upcase_column_fh(ofh, filenames):
    for fileno, filename in enumerate(filenames):
        with open(filename) as fh:
            for iline, line in enumerate(fh):
                if iline == 0:
                    if fileno == 0:
                        print(line, end='', file=ofh)
                    continue
                line = line.rstrip()
                F = line.split('\t')
                F[0] = F[0].upper()
                print(*F, sep='\t', file=ofh)


def gather_upcase_column(outfile, filenames):
    with open(outfile, 'w') as ofh:
        gather_upcase_column_fh(ofh, filenames)


def gather_and_split_fh(fh, o_microbiome, o_expression):
    data_type_names = DATA_TYPE_NAMES
    raw_data_types = ['kraken', 'htseq']
    ofhs = [o_microbiome, o_expression]
    lines = list(fh)
    for raw_data_type, ofh in zip(raw_data_types, ofhs):
        for lineno, line in enumerate(lines):
            if lineno == 0:
                print(line, end='', file=ofh)
                continue
            line = line.rstrip()
            F = line.split('\t')
            if F[2] != raw_data_type:
                continue
            F[2] = data_type_names.get(raw_data_type, F[2])
            print(*F, sep="\t", file=ofh)


def gather_and_split(filename, f_microbiome, f_expression):
    with (open(filename) as fh, open(f_microbiome, 'w') as
          o_microbiome, open(f_expression, 'w') as o_expression):
        gather_and_split_fh(fh, o_microbiome, o_expression)


def gather_s4_panel_fh(fh, filesi, filesii):
    gather_upcase_column_fh(fh, filesi)
    print('\n', file=fh)
    gather_upcase_column_fh(fh, filesii)


def gather_s4_panel(filename, filesi, filesii):
    with open(filename, 'w') as fh:
        gather_s4_panel_fh(fh, filesi, filesii)


if True:
    gather_separate_stats('figure2a.tsv', FIG2A)
    gather_separate_stats('figure2b.tsv', FIG2B)
    gather_upcase_column('figure2c.tsv', FIG2C)

    gather_separate_stats('figure3a.tsv', FIG3A)
    gather_upcase_column('figure3b.tsv', FIG3B)
    gather_upcase_column('figure3c.tsv', FIG3C)
    gather_parse_filename('figure3d.tsv', FIG3D)
    gather_parse_filename('figure3e.tsv', FIG3E)

    gather_separate_stats('figure4a.tsv', FIG4A)
    gather_upcase_column('figure4b.tsv', FIG4B)
    gather_upcase_column('figure4c.tsv', FIG4C)
    gather_parse_filename('figure4d.tsv', FIG4D)
    gather_parse_filename('figure4e.tsv', FIG4E)

    gather_and_split(FIG5AC, 'figure5a.tsv', 'figure5c.tsv')
    gather_and_split(FIG5BD, 'figure5b.tsv', 'figure5d.tsv')

    gather_parse_filename('figure6a.tsv', FIG6A)
    gather_parse_filename('figure6b.tsv', FIG6B)
    gather_parse_filename('figure6c.tsv', FIG6C)
    gather_parse_filename('figure6d.tsv', FIG6D)

    gather_upcase_column('figureS1a.tsv', FIGS1A)
    gather_parse_filename('figureS1b.tsv', FIGS1B)

    gather_upcase_column('figureS2a.tsv', FIGS2A)
    gather_parse_filename('figureS2b.tsv', FIGS2B)

    gather_upcase_column('figureS3a.tsv', FIGS3A)
    gather_parse_filename('figureS3b.tsv', FIGS3B)

    gather_s4_panel("figureS4a.tsv", FIGS4AI, FIGS4AII)
    gather_upcase_column('figureS4b.tsv', FIGS4B)
    gather_s4_panel('figureS4c.tsv', FIGS4CI, FIGS4CII)

    gather_parse_filename('figureS5a.tsv', FIGS5A)
    gather_parse_filename('figureS5b.tsv', FIGS5B)

    gather_parse_filename('figureS6a.tsv', FIGS6A)
    gather_parse_filename('figureS6b.tsv', FIGS6B)
