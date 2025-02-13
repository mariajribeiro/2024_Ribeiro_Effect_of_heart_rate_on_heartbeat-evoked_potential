% for each participant test the effect of task - Level 1 LIMO
% analysis of ECG amplitude
clear; close all;
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

older = [11 12 14 20 21 22 23 32 37 38 41 43 47 48 49 52 55 57 58 63 64 65 67 69 7 70 71 75 8 83 86];
young = [13 15 16 25 26 28 31 33 34 36 4 42 44 45 46 50 51 53 54 56 59 6 62 66 68 72 74 76 78 80 82 84 85 9];
% all files
files_folder = 'D:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\eeg_derivatives_05Hzfilter';
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
ECG4plot_avg = cell(2, 2); % group x task
yng = 0; old = 0;
for s = 1:length(subjects)

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
    
    ECG_data = [];
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
        
        % make epoch with all R-peaks
        EEG = pop_epoch( EEG, {  'Rpeak'  }, [-0.15         1.5], 'newname', 'Rpeak epochs all', 'epochinfo', 'yes');
%         EEG = pop_rmbase( EEG, [-150 -50] ,[]);
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
        
        % calculate heart rate for each epoch
        IBI = [];
        for ep = 1:length(EEG.epoch)
            Rpeaks = find(strcmp(EEG.epoch(ep).eventtype, 'Rpeak'));
            if length(Rpeaks) > 1
                IBI(ep) = EEG.epoch(ep).eventlatency{Rpeaks(2)};
            else
                IBI(ep) = NaN;
            end
        end
                
        % find outliers - ectopic beats
        zscore_IBI = (IBI - mean(IBI, 'omitnan'))/std(IBI, 'omitnan');
        IBI(abs(zscore_IBI) > 4) = NaN;
        
        % ECG excluding outliers that might refer to ectopic beats or
        % trials where the Rpeak was wrongly detected
        ECG_data = cat(3, ECG_data, EEG.data(60, 125:425, ~isnan(IBI)));
        number_trials(t) = size(EEG.data(:, :, ~isnan(IBI)), 3);
        
        if group == 1
            ECG4plot_avg{group, t}(yng, :) = squeeze(mean(EEG.data(60, 1:425, ~isnan(IBI)), 3));
        else
            ECG4plot_avg{group, t}(old, :) = squeeze(mean(EEG.data(60, 1:425, ~isnan(IBI)), 3));
        end
        
    end
    
    % create variable with data electrodes X time frames X trials -
    % starting at 
    % categorical variable with different task conditions
    categorical_variable=[ones(number_trials(1), 1); ones(number_trials(2), 1)*2];

    LIMO_directory=strcat(pwd, '\AB', num2str(subj_number));
    mkdir(LIMO_directory);
    
    save(strcat(LIMO_directory, filesep, 'ECG_data.mat'), 'ECG_data');
    save(strcat(LIMO_directory, filesep, 'categorical_variable.mat'), 'categorical_variable');
    
    % create LIMO variable
    LIMO.Level                    = 1;
    
    LIMO.Type = 'Channels';
    LIMO.Analysis = 'Time';
    
    LIMO.dir                      = LIMO_directory;

    LIMO.data.data_dir            = LIMO_directory;
    LIMO.data.data                = 'ECG_data.mat';

    LIMO.data.start               = 100;
    LIMO.data.end                 = 700;
    LIMO.data.trim1               = 1; 
    LIMO.data.trim2               = size(ECG_data, 2);
    
    LIMO.data.timevect            = 100:2:700;
    
    LIMO.data.sampling_rate       = 500;
    
    LIMO.data.Cat                 = categorical_variable; 
    LIMO.data.Cont                = [];

    load('D:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\LIMO_stats\expected_chanlocs_both.mat');

    LIMO.data.neighbouring_matrix = channeighbstructmat(1, 1);
    LIMO.data.chanlocs = expected_chanlocs(1);

    LIMO.design.fullfactorial     = 0;
    LIMO.design.zscore            = 0; % IT LOOKS LIKE IT IS ZSCORING THE DATA - CHECK!!!
    LIMO.design.method            = 'OLS'; 
    LIMO.design.type_of_analysis  = 'Mass-univariate';
    LIMO.design.bootstrap         = 0; 
    LIMO.design.tfce              = 0;

    [LIMO.design.X, LIMO.design.nb_conditions, LIMO.design.nb_interactions, LIMO.design.nb_continuous] = limo_design_matrix(ECG_data, LIMO, 0);

    LIMO.design.status          = 'to do';
    % LIMO.design.name  = 'paired t-test all electrodes';

    save(strcat(LIMO_directory, filesep, 'LIMO.mat'), 'LIMO');

    % run LIMO first level
    cd(LIMO_directory);
    limo_eeg(4)
    
    cd ..

end


save ECG4plot_avg ECG4plot_avg


