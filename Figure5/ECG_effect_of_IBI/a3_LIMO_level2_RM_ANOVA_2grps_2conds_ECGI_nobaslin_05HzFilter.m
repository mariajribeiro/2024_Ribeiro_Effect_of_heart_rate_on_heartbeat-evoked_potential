% run second level LIMO stats

clear;
eeglab; limo_eeg
% check if betas from regression of IBI vs HEP are significantly different
% from zero for passive, simple RT and go/no-go separately
% then check if betas are different across tasks
% then check if betas are different across groups
older = [11 12 14 20 21 22 23 32 37 38 41 43 47 48 49 52 55 57 58 63 64 65 67 69 7 70 71 75 8 83 86];
young = [13 15 16 25 26 28 31 33 34 36 4 42 44 45 46 50 51 53 54 56 59 6 62 66 68 72 74 76 78 80 82 84 85 9];
tasks = {'simpleRT', 'gng'};
folder_list = dir(pwd);
load('D:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\LIMO_stats\expected_chanlocs_both.mat');

level2_RM_ANOVA_dir = [pwd, '\level2_RM_ANOVA_IBIvsECG_50-350ms'];
mkdir(level2_RM_ANOVA_dir);
% for regression with IBI/HeartRate simpleRT, gng - beta = 4:6

cd(level2_RM_ANOVA_dir);
% limo_random_robust(6,y,gp,factor_levels,LIMO,'go',option)
%                    6 = Repeated measures ANOVA/ANCOVA using multivariate approach
%                    y = data (dim channels, time or freq, subjects, measures)
%                      = data (dim channels, freq, time, subjects, measures)
%                      = the name of the Yr file
%                    gp = a vector defining gps
%                    factor_levels = a vector specifying the levels of each repeated measure factor
%                    LIMO the basic structure with data, design and channel info
%                         or the full name of the LIMO file
%                    'go' is optional and prompt or not the pseudo-design
%                         options are 'yes' (prompt, default),
%                         or 'no' (no prompt, usuful for scripting)

s = 0; Y = [];
for f = 1:length(folder_list)
    if contains(folder_list(f).name, 'AB')
        s = s + 1;
        load([folder_list(f).folder, filesep, folder_list(f).name, filesep, 'Betas.mat'])
        Y(:, :, s, :) = Betas(:, 1:150, 3:4); %channels, time, subjects, measures
        subjects{s} = [folder_list(f).name];
        if ismember(str2num(folder_list(f).name(3:end)), young)
            gp(s, 1) = 1;
        elseif ismember(str2num(folder_list(f).name(3:end)), older)
            gp(s, 1) = 2;
        end
    end
end

save Y Y
save gp gp

% create LIMO variable
LIMO.Level                    = 2;

LIMO.Type = 'Channels';
LIMO.Analysis = 'Time';

LIMO.dir                      = level2_RM_ANOVA_dir;

LIMO.data.data_dir            = level2_RM_ANOVA_dir;
% LIMO.data.data                = 

LIMO.data.start               = 51;
LIMO.data.end                 = 400;
LIMO.data.trim1               = 1; 
LIMO.data.trim2               = size(Y, 2);
LIMO.data.timevect            = 51:2:400;
LIMO.data.sampling_rate       = 500;
LIMO.data.neighbouring_matrix = channeighbstructmat(1, 1);
LIMO.data.chanlocs = expected_chanlocs(1);

% LIMO.design.fullfactorial     = 0;
% LIMO.design.zscore            = 0;
LIMO.design.method            = 'Trimmed mean'; 
LIMO.design.type_of_analysis  = 'Mass-univariate';
LIMO.design.bootstrap         = 1000; 
LIMO.design.tfce              = 0;

save LIMO LIMO

%% repeated measures ANOVA
% limo_random_robust(6,y,gp,factor_levels,LIMO,'go',option)
LIMOpath = limo_random_robust(6, Y, gp, 2, LIMO);
% 1 one_sample_parameter_X (channels, frames [time, freq or freq-time], [mean value, se, df, t, p])
%   H0_one_sample_ttest_parameter_X (channels, frames, [T values under H0, p values under H0], LIMO.design.bootstrap)

%% view results  - clustering
% addpath('C:\eeglab2022.0\plugins\limo_tools-master_Nov2022') 
% level2_RM_ANOVA_dir = [pwd, '\level2_RM_ANOVA_IBIvsHEP'];
% cd(level2_RM_ANOVA_dir);
% https://github.com/LIMO-EEG-Toolbox/limo_tools/wiki/Workspace_variables
mask = [];
LIMO_file = fullfile(pwd,'LIMO.mat');
% limo_display_results(Type,FileName,PathName,p,MCC,LIMO,flag,options)
limo_display_results(1,'Rep_ANOVA_Main_effect_1.mat',pwd,0.05, 2,...
    LIMO_file, 0);

if ~isempty(mask)
    save mask_main_effect mask
    save stat_values_main_effect stat_values %F/t values
    save p_values_main_effect p_values


    number_clusters = max(unique(mask(:)));
    p_value_cluster = p_values(mask == 1);
    max_p_value = max(p_values(mask ~= 0));
    %     what is the range of significant F/t values? 
    stats_F = [min(stat_values(mask ~= 0)) max(stat_values(mask ~= 0))];
end

%%
limo_display_results(1,'Rep_ANOVA_Interaction_gp_Factor_1.mat',pwd,0.05, 2,...
    fullfile(pwd,'LIMO.mat'), 0);

limo_display_results(1,'Rep_ANOVA_Gp_effect.mat',pwd,0.05, 2,...
    fullfile(pwd,'LIMO.mat'), 0);
%% saveas(gcf, 'One_sample_timecourse.fig'); close(gcf)
if ~isempty(mask)
    save mask mask
    save stat_values stat_values %F/t values
    save p_values p_values
end

%     number_clusters = max(unique(mask(:)));
%     p_value_cluster = p_values(mask == 1);
%     max_p_value = max(p_values(mask ~= 0));
%     %     what is the range of significant F/t values? 
%     % [min(stat_values(mask ~= 0)) max(stat_values(mask ~= 0))] etc ..
% 
%     % replot using 
%     limo_display_image(LIMO,stat_values ,mask,'go/no-go HR effect',0)  


%% effect sizes
level2_RM_ANOVA_dir = [pwd, '\level2_RM_ANOVA_IBIvsECG_50-350ms'];
cd(level2_RM_ANOVA_dir);

clear

load LIMO.mat
LIMO.design.name = 'Repeated Measures ANOVA';
LIMO.design.effects = 'Interaction_gp_Factor_1';
save LIMO LIMO

% group effect
[name] = limo_get_effect_size('Rep_ANOVA_Gp_effect.mat')
load('Rep_ANOVA_Gp_effect_PartialEta2.mat')
data_effect_size = reshape(effect_size, size(effect_size, 1)*size(effect_size, 2), 1);
max(data_effect_size)
median(data_effect_size)

clear; clc
% interaction effect
[name] = limo_get_effect_size('Rep_ANOVA_Interaction_gp_Factor_1.mat')
load(name)
data_effect_size = reshape(effect_size, size(effect_size, 1)*size(effect_size, 2), 1);
max(data_effect_size)
median(data_effect_size)

% task effect
[name] = limo_get_effect_size('Rep_ANOVA_Main_effect_1.mat')
load(name)
data_effect_size = reshape(effect_size, size(effect_size, 1)*size(effect_size, 2), 1);
max(data_effect_size)
median(data_effect_size)

