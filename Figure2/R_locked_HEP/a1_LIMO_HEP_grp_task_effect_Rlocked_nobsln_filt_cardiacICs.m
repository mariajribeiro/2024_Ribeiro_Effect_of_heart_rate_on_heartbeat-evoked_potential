% study HEP locked with the R peak R-HEP
% where it is significantly different from zero
% where there are group effects
clear; close all;
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

addpath('E:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study');

older = [11 12 14 20 21 22 23 32 37 38 41 43 47 48 49 52 55 57 58 63 64 65 67 69 7 70 71 75 8 83 86];
young = [13 15 16 25 26 28 31 33 34 36 4 42 44 45 46 50 51 53 54 56 59 6 62 66 68 72 74 76 78 80 82 84 85 9];
% all files
files_folder = 'E:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\eeg_derivatives_05Hzfilt_cardiacICs';
folder_content = dir(files_folder);

% subjects included
subjects = {};
for f = 1:length(folder_content)
    if contains(folder_content(f).name, 'eeg_Rpeak_Tpeak')
        subjects = [subjects, folder_content(f).name(21:24)];
    end
end
subjects = unique(subjects);

tasks = {'simpleRT', 'gonogo'};%tasks = {'passive', 'simpleRT', 'gonogo'};
HEP4plot_avg = cell(2, 2); % group x task
yng = 0; old = 0;
for s =  1:length(subjects)

    if strcmp(subjects{s}(4), '_')
        sbj = 3;
    else
        sbj = 3:4;
    end
    
    subj_number = str2num(subjects{s}(sbj));
    if ismember(subj_number, young)
        group = 1; yng = yng + 1;
    else
        group = 2; old = old + 1;
    end
    
    HEP_data = [];
    IBI = cell(3, 1); number_trials = [];
       
    for t = 1:length(tasks)
        % clear eeglab
        STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];  

        % load data
        filename = {};
        for f = 1:length(folder_content)
            if contains(folder_content(f).name, subjects(s)) && contains(folder_content(f).name, tasks(t))
                filename = folder_content(f).name;
                EEG = pop_loadset('filename', filename, 'filepath', files_folder );
%                 EEG = pop_select( EEG, 'nochannel',{'EKG'});
                [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
            end
        end

        % append data from both runs
        if length(ALLEEG) == 2
            EEG = pop_mergeset( ALLEEG, [1  2], 0);
            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'gui','off'); 
        end
        
        %% delete R and T peaks associated with abnormal cardiac cycles or
        % cycles where the cycle was not segmented properly
        outliers_evnt = find_cardiac_cycle_outliers(EEG);
        %% delete events associated with outlier cycles
        EEG = pop_selectevent( EEG, 'omitevent', outliers_evnt ,'deleteevents','on');
        
        % make epoch with all R-peaks
        EEG = pop_epoch( EEG, {  'Rpeak'  }, [-0.15   .7], 'newname', 'Rpeak epochs all', 'epochinfo', 'yes');
%         EEG = pop_rmbase( EEG, [-150 -50] ,[]);
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
        % analyse data from 300 ms up to 650 ms after R peak
        HEP_data = cat(3, HEP_data, EEG.data(1:59, 226:400, :));
        number_trials(t) = size(EEG.data, 3);
        
        if group == 1
            HEP4plot_avg{group, t}(yng, :, :) = mean(EEG.data, 3);
        else
            HEP4plot_avg{group, t}(old, :, :) = mean(EEG.data, 3);
        end
        
    end

    % create variable with data electrodes X time frames X trials -
    % starting at 
    % categorical variable with different task conditions
    categorical_variable=[ones(number_trials(1), 1); ones(number_trials(2), 1)*2];

    LIMO_directory=strcat(pwd, '\AB', num2str(subj_number));
    mkdir(LIMO_directory);
    
    save(strcat(LIMO_directory, filesep, 'HEP_data.mat'), 'HEP_data');
    save(strcat(LIMO_directory, filesep, 'categorical_variable.mat'), 'categorical_variable');
    
    % create LIMO variable
    LIMO.Level                    = 1;
    
    LIMO.Type = 'Channels';
    LIMO.Analysis = 'Time';
    
    LIMO.dir                      = LIMO_directory;

    LIMO.data.data_dir            = LIMO_directory;
    LIMO.data.data                = 'HEP_data.mat';
    % analyse data from 300 ms up to 650 ms after R peak
    LIMO.data.start               = 301;
    LIMO.data.end                 = 649;
    LIMO.data.trim1               = 1; 
    LIMO.data.trim2               = size(HEP_data, 2);
    
    LIMO.data.timevect            = 301:2:649;
    
    LIMO.data.sampling_rate       = 500;
    
    LIMO.data.Cat                 = categorical_variable; 
    LIMO.data.Cont                = [];

    load('E:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\LIMO_stats\expected_chanlocs_both.mat');

    LIMO.data.neighbouring_matrix = channeighbstructmat;
    LIMO.data.chanlocs = expected_chanlocs;

    LIMO.design.fullfactorial     = 0;
    LIMO.design.zscore            = 0; % IT LOOKS LIKE IT IS ZSCORING THE DATA - CHECK!!!
    LIMO.design.method            = 'OLS'; 
    LIMO.design.type_of_analysis  = 'Mass-univariate';
    LIMO.design.bootstrap         = 0; 
    LIMO.design.tfce              = 0;

    [LIMO.design.X, LIMO.design.nb_conditions, LIMO.design.nb_interactions, LIMO.design.nb_continuous] = limo_design_matrix(HEP_data, LIMO, 0);

    LIMO.design.status          = 'to do';
    % LIMO.design.name  = 'paired t-test all electrodes';

    save(strcat(LIMO_directory, filesep, 'LIMO.mat'), 'LIMO');

    % run LIMO first level
    cd(LIMO_directory);
    limo_eeg(4)
    
    cd ..

end

save HEP4plot_avg HEP4plot_avg
