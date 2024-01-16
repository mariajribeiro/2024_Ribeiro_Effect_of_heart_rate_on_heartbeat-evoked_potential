% LIMO level 1 effect of task and group on T-locked HEP - with random
% events
clear; close all;
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

older = [11 12 14 20 21 22 23 32 37 38 41 47 48 49 52 55 57 58 63 64 65 67 69 7 70 71 75 8 83 86];
young = [13 15 16 25 26 28 31 33 34 36 4 42 44 45 46 50 51 53 54 56 59 6 62 66 68 72 74 76 78 80 82 84 85 9];
% note: AB21, AB50 and AB55 were excluded from T peak analyses because (s = 8 30 35)
% these were not clearly visible in the ECGs

% load RT intervals for each participant
RT = readtable('D:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\LIMO_stats\HEP_grp_task_effect_Tlocked\RT_interval.xlsx');
   
% all files
files_folder = 'D:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\eeg_derivatives\random_R_T_events';
folder_content = dir(files_folder);

% subjects included
subjects = {};
for f = 1:length(folder_content)
    if contains(folder_content(f).name, 'eeg_random_Rpeak_Tpeak')
        subjects = [subjects, folder_content(f).name(28:31)];
    end
end
subjects = unique(subjects);
% subjects excluded from T-locked HEP analyses due to difficulty
% identifying AB21, AB50 and AB55
subjects = setdiff(subjects, {'AB21', 'AB50', 'AB55'});

tasks = {'simpleRT', 'gonogo'};

HEP_T_locked_random = cell(2, 2);

number_trials_random_T = cell(2); % group - number of trials: all, after cue, b4 cue
number_trials_random_T{1} = zeros(length(young), 3);
number_trials_random_T{2} = zeros(length(older), 3);

yng = 0; old = 0; without_t = [];
for s = 1:length(subjects)
    if strcmp(subjects{s}(4), '_')
        sbj = 3;
    else
        sbj = 3:4;
    end
    
    subj_number = str2double(subjects{s}(sbj));
    if ismember(subj_number, young)
        group = 1; yng = yng + 1;
    else
        group = 2; old = old + 1;
    end
    
    HEP_T_locked_random_data = [];
    number_trials = [];
       
    for t = 1:length(tasks)
        % clear eeglab
        STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];  

        % load data
        filename = {};
        for f = 1:length(folder_content)
            if contains(folder_content(f).name, subjects(s)) && contains(folder_content(f).name, tasks(t))
                filename = folder_content(f).name;
                EEG = pop_loadset('filename', filename, 'filepath', files_folder );
                EEG = pop_select( EEG, 'nochannel',{'EKG'}); % INCLUDE EKG IN ANALYSES!!!!
                [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
            end
        end

        % append data from both runs
        if length(ALLEEG) == 2
            EEG = pop_mergeset( ALLEEG, [1  2], 0);
            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'gui','off'); 
        end
        
        % make epoch with all R-peaks
        EEG = pop_epoch( EEG, {  'Rpeak'  }, [-0.15         1.5], 'newname', 'Rpeak epochs all', 'epochinfo', 'yes');
        EEG = pop_rmbase( EEG, [-150 -50] ,[]);
        RT_interval = [];
        if t == 1
            RT_interval = mean([RT.RT_interv_simpleRT(RT.participant == subj_number), RT.RT_interv_gng(RT.participant == subj_number)]);
        elseif t == 2
            RT_interval = RT.RT_interv_simpleRT(RT.participant == subj_number);
        elseif t == 3
            RT_interval = RT.RT_interv_gng(RT.participant == subj_number);
        end
        % cut epoch locked to R random plus the avg RT interval of that
        % participant
        EEG = pop_select( EEG, 'time',[(RT_interval+150-200)*.001  (RT_interval+150+400)*.001]);
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );

        % HEP T-locked data
        HEP_T_locked_random_data = cat(3, HEP_T_locked_random_data, EEG.data(1:59, 126:end, :));
        number_trials(t) = size(EEG.data, 3);
        
        if group == 1
            HEP_T_locked_random{group, t}(yng, :, :) = mean(EEG.data, 3);
        else
            HEP_T_locked_random{group, t}(old, :, :) = mean(EEG.data, 3);
        end
    end

    if sum(number_trials) > 0
        % create variable with data electrodes X time frames X trials -
        % starting at 
        % categorical variable with different task conditions
        categorical_variable=[ones(number_trials(1), 1); ones(number_trials(2), 1)*2];

        LIMO_directory=strcat(pwd, '\AB', num2str(subj_number));
        mkdir(LIMO_directory);

        save(strcat(LIMO_directory, filesep, 'HEP_T_locked_random_data.mat'), 'HEP_T_locked_random_data');
        save(strcat(LIMO_directory, filesep, 'categorical_variable.mat'), 'categorical_variable');

%         % create LIMO variable
%         LIMO.Level                    = 1;
% 
%         LIMO.Type = 'Channels';
%         LIMO.Analysis = 'Time';
% 
%         LIMO.dir                      = LIMO_directory;
% 
%         LIMO.data.data_dir            = LIMO_directory;
%         LIMO.data.data                = 'HEP_T_locked_random_data.mat';
% 
%         LIMO.data.start               = 51;
%         LIMO.data.end                 = 600;
%         LIMO.data.trim1               = 1; 
%         LIMO.data.trim2               = size(HEP_T_locked_random_data, 2);
% 
%         LIMO.data.timevect            = 51:2:600;
% 
%         LIMO.data.sampling_rate       = 500;
% 
%         LIMO.data.Cat                 = categorical_variable; 
%         LIMO.data.Cont                = [];
% 
%         load('D:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\LIMO_stats\expected_chanlocs_both.mat');
% 
%         LIMO.data.neighbouring_matrix = channeighbstructmat;
%         LIMO.data.chanlocs = expected_chanlocs;
% 
%         LIMO.design.fullfactorial     = 0;
%         LIMO.design.zscore            = 0; % IT LOOKS LIKE IT IS ZSCORING THE DATA - CHECK!!!
%         LIMO.design.method            = 'OLS'; 
%         LIMO.design.type_of_analysis  = 'Mass-univariate';
%         LIMO.design.bootstrap         = 0; 
%         LIMO.design.tfce              = 0;
% 
%         [LIMO.design.X, LIMO.design.nb_conditions, LIMO.design.nb_interactions, LIMO.design.nb_continuous] = limo_design_matrix(HEP_T_locked_random_data, LIMO, 0);
% 
%         LIMO.design.status          = 'to do';
%         % LIMO.design.name  = 'paired t-test all electrodes';
% 
%         save(strcat(LIMO_directory, filesep, 'LIMO.mat'), 'LIMO');
% 
%         % run LIMO first level
%         cd(LIMO_directory);
%         limo_eeg(4)
% 
%         cd ..
    end
end
save HEP_T_locked_random HEP_T_locked_random