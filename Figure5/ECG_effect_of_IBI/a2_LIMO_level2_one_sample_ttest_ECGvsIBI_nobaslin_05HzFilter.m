% run second level LIMO stats
% addpath LIMO_eeg
clear; close all
eeglab; limo_eeg;
% 
% check if betas from regression of IBI vs HEP are significantly different
% from zero for all trials (simple_RT + go/no-go)

tasks = {'simpleRT', 'gng'};
folder_list = dir(pwd);

load('D:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\LIMO_stats\expected_chanlocs_both.mat');
level2_onesample_dir = [pwd filesep 'level2_onesample_analysis_50-350ms'];
mkdir(level2_onesample_dir);
% for regression with IBI/HeartRate simpleRT, gng - beta = 3:4
for t = 1:2 % two tasks

    LIMO_directory = [level2_onesample_dir filesep tasks{t}];
    mkdir(LIMO_directory);
    cd(LIMO_directory);
    % create variable with data Y = channel x time x subject

    s = 0; Y = [];
    for f = 1:length(folder_list)
        if contains(folder_list(f).name, 'AB')
            s = s + 1;
            load([folder_list(f).folder, filesep, folder_list(f).name, filesep, 'Betas.mat'])
            Y(:, :, s) = squeeze(Betas(:, 1:150, t+2));
            subjects{s} = [folder_list(f).name];
        end
    end

    save Y Y

    % create LIMO variable
    LIMO.Level                    = 2;

    LIMO.Type = 'Channels';
    LIMO.Analysis = 'Time';

    LIMO.dir                      = LIMO_directory;

    LIMO.data.data_dir            = LIMO_directory;
    % LIMO.data.data                = 

    LIMO.data.start               = 51;
    LIMO.data.end                 = 349;
    LIMO.data.trim1               = 1; 
    LIMO.data.trim2               = size(Y, 2);
    LIMO.data.timevect            = 51:2:349;
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

    % one-sample t-test to check where correlation between heart rate and HEP
    % amplitude is significant

    LIMOpath = limo_random_robust(1, Y, 1, LIMO);
    % 1 one_sample_parameter_X (channels, frames [time, freq or freq-time], [mean value, se, df, t, p])
    %   H0_one_sample_ttest_parameter_X (channels, frames, [T values under H0, p values under H0], LIMO.design.bootstrap)

    %% view results  - clustering
    % https://github.com/LIMO-EEG-Toolbox/limo_tools/wiki/Workspace_variables
    mask = [];
    % limo_display_results(Type,FileName,PathName,p,MCC,LIMO,flag,options)
    limo_display_results(1,'one_sample_ttest_parameter_1.mat',pwd,0.05, 2,...
        fullfile(pwd,'LIMO.mat'), 0);
    % saveas(gcf, 'One_sample_timecourse.fig'); close(gcf)
    if ~isempty(mask)
        save mask mask
        save stat_values stat_values %F/t values
        save p_values p_values


        number_clusters = max(unique(mask(:))); 
        p_value_cluster = p_values(mask == 1);
        max_p_value = max(p_values(mask ~= 0));
        %     what is the range of significant F/t values? 
        stats_F = [min(stat_values(mask ~= 0)) max(stat_values(mask ~= 0))];
    end
%     % replot using 
%     limo_display_image(LIMO,stat_values ,mask,'go/no-go HR effect',0)  

end

%% effect sizes
tasks = {'simpleRT', 'gng'};
level2_onesample_dir = [pwd filesep 'level2_onesample_analysis_50-350ms'];
t = 2;
LIMO_directory = [level2_onesample_dir filesep tasks{t}];
cd(LIMO_directory);

load LIMO.mat
LIMO.data.neighbouring_matrix = 1;
save LIMO LIMO

clear
load mask
[name] = limo_get_effect_size('one_sample_ttest_parameter_1.mat')
load(name)
max(effect_size(mask > 0))
median(effect_size(mask > 0))

data_effect_size = reshape(effect_size, size(effect_size, 1)*size(effect_size, 2), 1);
max(data_effect_size)
median(data_effect_size)

