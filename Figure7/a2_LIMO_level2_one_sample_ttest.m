% run second level LIMO stats

clear; close all
% addpath LIMO_eeg
eeglab
limo_eeg;
% 
% check if betas from regression of IBI vs HEP are significantly different
% from zero for simple RT and go/no-go separately
% then check if betas are different across tasks
% then check if betas are different across groups

tasks = {'passive', 'simpleRT', 'gng'};
folder_list = dir(pwd);
load('D:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\LIMO_stats\expected_chanlocs_both.mat');

level2_onesample_dir = [pwd, '\level2_onesample_IBIvsHEP_3tasks'];
mkdir(level2_onesample_dir);
% for regression with IBI/HeartRate simpleRT, gng - beta = 3:4
for t = 1:length(tasks) % task conditions

    LIMO_directory = [level2_onesample_dir filesep tasks{t}];
    mkdir(LIMO_directory);
    cd(LIMO_directory);
    % create variable with data Y = channel x time x subject

    s = 0; Y = [];
    for f = 1:length(folder_list)
        if contains(folder_list(f).name, 'AB')
            s = s + 1;
            load([folder_list(f).folder, filesep, folder_list(f).name, filesep, 'Betas.mat'])
            Y(:, :, s) = squeeze(Betas(:, :, t+3));
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



