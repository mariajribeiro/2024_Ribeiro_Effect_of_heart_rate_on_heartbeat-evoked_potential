% run second level LIMO stats
% addpath LIMO_eeg
clear; close all

% check where the HEP is significantly different from zero - all
older = [11 12 14 20 21 22 23 32 37 38 41 43 47 48 49 52 55 57 58 63 64 65 67 69 7 70 71 75 8 83 86];
young = [13 15 16 25 26 28 31 33 34 36 4 42 44 45 46 50 51 53 54 56 59 6 62 66 68 72 74 76 78 80 82 84 85 9];
groups = {young, older};
grp_name = {'young', 'older'};
folder_list = dir(pwd);
load('E:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\LIMO_stats\expected_chanlocs_both.mat');
group = {'young', 'older'};
level2_onesample_dir = [pwd, '\level2_onesample_HEP_grps_together'];
mkdir(level2_onesample_dir);

% create variable with data Y = channel x time x subject

s = 0; Y = []; subjects = {};
for f = 1:length(folder_list)
    if contains(folder_list(f).name, 'AB')
        if ismember(str2double(folder_list(f).name(3:end)), groups{1})
            load([folder_list(f).folder, filesep, folder_list(f).name, filesep, 'HEP_T_locked_data.mat'])
            s = s + 1;
            Y(:, :, s) = squeeze(mean(HEP_T_locked_data(:, 126:300, :), 3)); % channels x time x subjects - from 51 to 400 ms locked with T peak
            subjects{s} = [folder_list(f).name];
        elseif ismember(str2double(folder_list(f).name(3:end)), groups{2})
            load([folder_list(f).folder, filesep, folder_list(f).name, filesep, 'HEP_T_locked_data.mat'])
            s = s + 1;
            Y(:, :, s) = squeeze(mean(HEP_T_locked_data(:, 126:300, :), 3)); % channels x time x subjects - from 51 to 400 ms locked with T peak
            subjects{s} = [folder_list(f).name];
        end
    end
end
%%
cd(level2_onesample_dir)
% create LIMO variable
LIMO.Level                    = 2;

LIMO.Type = 'Channels';
LIMO.Analysis = 'Time';

LIMO.dir                      = level2_onesample_dir;

LIMO.data.data_dir            = level2_onesample_dir;
% LIMO.data.data                = 

LIMO.data.start               = 51;
LIMO.data.end                 = 400;
LIMO.data.trim1               = 1; 
LIMO.data.trim2               = size(Y, 2);
LIMO.data.timevect            = 51:2:400;
LIMO.data.sampling_rate       = 500;
LIMO.data.neighbouring_matrix = channeighbstructmat;
LIMO.data.chanlocs = expected_chanlocs;

% LIMO.design.fullfactorial     = 0;
% LIMO.design.zscore            = 0;
LIMO.design.method            = 'Trimmed mean'; 
LIMO.design.type_of_analysis  = 'Mass-univariate';
LIMO.design.bootstrap         = 1000; 
LIMO.design.tfce              = 0;

save LIMO LIMO
eeglab; limo_eeg;
% one-sample t-test to check where correlation between heart rate and HEP
% amplitude is significant
%%
LIMOpath = limo_random_robust(1, Y, 1, LIMO);
% 1 one_sample_parameter_X (channels, frames [time, freq or freq-time], [mean value, se, df, t, p])
%   H0_one_sample_ttest_parameter_X (channels, frames, [T values under H0, p values under H0], LIMO.design.bootstrap)

%% view results  - clustering
% https://github.com/LIMO-EEG-Toolbox/limo_tools/wiki/Workspace_variables

% color maps from http://www.fabiocrameri.ch/colourmaps.php

mask = [];
% limo_display_results(Type,FileName,PathName,p,MCC,LIMO,flag,options)
limo_display_results(1,'one_sample_ttest_parameter_1.mat',pwd,0.05, 2,...
    fullfile(pwd,'LIMO.mat'), 0);
% saveas(gcf, 'One_sample_timecourse.fig'); close(gcf)
if ~isempty(mask)
    save mask mask
    save stat_values stat_values %F/t values
    save p_values p_values

    [row,col] = find(mask == 1);
    max_time = max(col)*2+100;

    number_clusters = max(unique(mask(:)));
    p_value_cluster = p_values(mask == 1);
    max_p_value = max(p_values(mask ~= 0));
    %     what is the range of significant F/t values? 
    stats_F = [min(stat_values(mask ~= 0)) max(stat_values(mask ~= 0))];

    % mean value in cluster
    load one_sample_ttest_parameter_1
    mean_data = squeeze(one_sample(:, :, 1));
    mean_value = mean(mean(mean_data(mask == 1)));

    % degress of freedom
    df = squeeze(one_sample(:, :, 3));
    t_value = squeeze(one_sample(:, :, 4));
end

%% replot using 
%     limo_display_image(LIMO,stat_values ,mask,'',0)  




