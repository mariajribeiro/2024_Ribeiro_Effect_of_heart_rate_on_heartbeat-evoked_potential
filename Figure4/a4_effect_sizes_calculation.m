% calculate effect sizes
eeglab
% [name,clusters] = limo_get_effect_size(file,mask)

cd([pwd filesep 'level2_RM_ANOVA_HEP_TLocked_b4_after_cue'])

load LIMO.mat

LIMO.design.name = 'Repeated Measures ANOVA';

save LIMO LIMO

% group effect
load mask_group_effect.mat
[name, cluster_groupeffectsize] = limo_get_effect_size('Rep_ANOVA_Gp_effect.mat', mask)

% main effect 1 - task
[name] = limo_get_effect_size('Rep_ANOVA_Main_effect_1.mat')


% main effect 2 - before/after
[name] = limo_get_effect_size('Rep_ANOVA_Main_effect_2.mat')
load('Rep_ANOVA_Main_effect_2_MahalanobisD.mat')
data_effect_size = reshape(effect_size, size(effect_size, 1)*size(effect_size, 2), 1);
max(data_effect_size)
median(data_effect_size)



