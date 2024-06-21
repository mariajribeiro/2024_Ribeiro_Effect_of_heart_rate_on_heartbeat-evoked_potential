% run second level LIMO stats
% test if diff in regression coefficients across tasks is related to
% differences in pupil phasic responses across tasks
clear;
eeglab; limo_eeg;

% addpath('C:\eeglab2022.0\plugins\limo_tools-master_Nov2022') 
% check if betas from regression of IBI vs HEP are correlated with
% amplitude of phasic pupil responses
older = [11 12 14 20 21 22 23 32 37 38 41 43 47 48 49 52 55 57 58 63 64 65 67 69 7 70 71 75 8 83 86];
young = [13 15 16 25 26 28 31 33 34 36 4 42 44 45 46 50 51 53 54 56 59 6 62 66 68 72 74 76 78 80 82 84 85 9];
tasks = {'simpleRT', 'gng'};
folder_list = dir(pwd);
load('D:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\LIMO_stats\expected_chanlocs_both.mat');

% load pupil data
pupil_dir = 'D:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\LIMO_stats\HEP_TLockedvsIBI';
load([pupil_dir, '\pupil_dilation_amplitude\BlnDivided_Part_PeakAmpLat_Median_D_CueLocked_OlderGrp.mat']);
load([pupil_dir, '\pupil_dilation_amplitude\BlnDivided_Part_PeakAmpLat_Median_G_CueLocked_OlderGrp']);
load([pupil_dir, '\pupil_dilation_amplitude\BlnDivided_Part_PeakAmpLat_Median_D_CueLocked_YoungGrp']);
load([pupil_dir, '\pupil_dilation_amplitude\BlnDivided_Part_PeakAmpLat_Median_G_CueLocked_YoungGrp']);

pupil_amplitude{1} = [BlnDivided_Part_PeakAmpLat_Median_D_CueLocked_YoungGrp(:, 1:2); ...
    BlnDivided_Part_PeakAmpLat_Median_D_CueLocked_OlderGrp(:, 1:2)];

pupil_amplitude{2} = [BlnDivided_Part_PeakAmpLat_Median_G_CueLocked_YoungGrp(:, 1:2); ...
    BlnDivided_Part_PeakAmpLat_Median_G_CueLocked_OlderGrp(:, 1:2)];

pupil_amplitude_diff = pupil_amplitude{2}(:, 2)-pupil_amplitude{1}(:, 2);



% limo_random_robust(4,y,X,parameter number,LIMO)
%                    4 = regression analysis
%                    y = data (dim channels, time or freq, subjects)
%                      = data (dim channels, freq, time, subjects)
%                    X = continuous regressor(s)
%                    parameter number = describe which parameter is currently analysed (e.g. 1 - use for maming only)



level2_regression_pupil_dir = [pwd, '\level2_regression_pupil_diff\'];
mkdir(level2_regression_pupil_dir);
% for regression with IBI/HeartRate simpleRT, gng - beta = 3:4

cd(level2_regression_pupil_dir);

s = 0; Y = []; pupil = [];
for f = 1:length(folder_list)
    if contains(folder_list(f).name, 'AB')
        % check if participant has pupil data
        subj_idx = find(pupil_amplitude{1}(:, 1) == str2num(folder_list(f).name(3:end)));
        if ~isempty(subj_idx)
            s = s + 1;
            pupil(s, 1) = pupil_amplitude_diff(subj_idx, 1);
            load([folder_list(f).folder, filesep, folder_list(f).name, filesep, 'Betas.mat'])
            % diff in regression coefficients across tasks
            Y(:, :, s) = Betas(:, 1:175, 4)-Betas(:, 1:175, 3); %channels, time, subjects)
            subjects{s} = [folder_list(f).name];
        end

    end
end

save Y Y

% create LIMO variable
LIMO.Level                    = 2;

LIMO.Type = 'Channels';
LIMO.Analysis = 'Time';

LIMO.dir                      = level2_regression_pupil_dir;

LIMO.data.data_dir            = level2_regression_pupil_dir;
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
%     LIMO.design.method            = 'Trimmed mean'; 
LIMO.design.type_of_analysis  = 'Mass-univariate';
LIMO.design.bootstrap         = 1000;
LIMO.design.tfce              = 0;

save LIMO LIMO

% regression

% limo_random_robust(4,y,X,parameter number,LIMO)
LIMOpath = limo_random_robust(4, Y, pupil, 1,LIMO);


%% view results  - clustering
% addpath('C:\eeglab2022.0\plugins\limo_tools-master_Nov2022') 
% level2_regression_pupil_dir = [pwd, '\level2_RM_ANOVA_IBIvsHEP'];
% cd(level2_regression_pupil_dir);
% https://github.com/LIMO-EEG-Toolbox/limo_tools/wiki/Workspace_variables
mask = [];
% limo_display_results(Type,FileName,PathName,p,MCC,LIMO,flag,options)
limo_display_results(1,'Covariate_effect_1.mat',pwd,0.05, 2,...
    fullfile(pwd,'LIMO.mat'), 0);

if ~isempty(mask)
    save mask_covariate mask
    save stat_covariate stat_values %F/t values
    save p_values_covariate p_values


    number_clusters = max(unique(mask(:)));
    p_value_cluster = p_values(mask == 1);
    max_p_value = max(p_values(mask ~= 0));
    %     what is the range of significant F/t values? 
    stats_F = [min(stat_values(mask ~= 0)) max(stat_values(mask ~= 0))];
end


