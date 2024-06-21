% run second level LIMO stats
% addpath LIMO_eeg
clear; close all

% check where the HEP is significantly different from zero - all
% participants together

folder_list = dir(pwd);
load('M:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\LIMO_stats\expected_chanlocs_both.mat');
group = {'young', 'older'};
level2_onesample_dir = [pwd, '\level2_onesample_HEP_all_partcts_50-400ms'];
mkdir(level2_onesample_dir);
% for regression with IBI/HeartRate passive, simpleRT, gng - beta = 4:6

LIMO_directory = level2_onesample_dir;
mkdir(LIMO_directory);
cd(LIMO_directory);
% create variable with data Y = channel x time x subject

s = 0; Y = []; subjects = {};
for f = 1:length(folder_list)
    if contains(folder_list(f).name, 'AB')
            load([folder_list(f).folder, filesep, folder_list(f).name, filesep, 'HEP_T_locked_random_data.mat'])
            s = s + 1;
            Y(:, :, s) = squeeze(mean(HEP_T_locked_random_data, 3)); % channels x time x subjects - from 51 to 400 ms locked with T peak
            subjects{s} = [folder_list(f).name];
    end
end

% create LIMO variable
LIMO.Level                    = 2;

LIMO.Type = 'Channels';
LIMO.Analysis = 'Time';

LIMO.dir                      = LIMO_directory;

LIMO.data.data_dir            = LIMO_directory;
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
% one-sample t-test to check where HEP
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
    limo_display_image(LIMO,stat_values ,mask,'',0)  




