% for each participant run a correlation between ECG and IBI on a
% trial-by-trial basis - including effect of task - only including simple
% RT and gng
clear; close all;
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
load('D:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\LIMO_stats\expected_chanlocs_both.mat');

older = [11 12 14 20 21 22 23 32 37 38 41 43 47 48 49 52 55 57 58 63 64 65 67 69 7 70 71 75 8 83 86];
young = [13 15 16 25 26 28 31 33 34 36 4 42 44 45 46 50 51 53 54 56 59 6 62 66 68 72 74 76 78 80 82 84 85 9];
% all files
files_folder = 'D:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\eeg_derivatives_05Hzfilter';
folder_content = dir(files_folder);

% ECG was inverted in some participants - need to invert before analyses
ECG_invert = [6 51 74 75 78];

% subjects included
subjects = {};
for f = 1:length(folder_content)
    if contains(folder_content(f).name, 'eeg_Rpeak_Tpeak')
        subjects = [subjects, folder_content(f).name(21:24)];
    end
end
subjects = unique(subjects);
% subjects excluded from T-locked HEP analyses due to difficulty
% identifying AB21, AB50 and AB55 and AB67
subjects = setdiff(subjects, {'AB21', 'AB50', 'AB55'});

tasks = {'simpleRT', 'gonogo'};

yng = 0; old = 0;
for s = 1:length(subjects)

    if strcmp(subjects{s}(4), '_')
        sbj = 3;
    else
        sbj = 3:4;
    end
    
    subj_number = str2num(subjects{s}(sbj));
%     if ismember(subj_number, young)
%         group = 1; yng = yng + 1;
%     else
%         group = 2; old = old + 1;
%     end
%     
    ECG_T_locked_data = []; IBI = cell(2, 1); number_trials = [];
       
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
                if ismember(subj_number, ECG_invert)
                    EEG.data(60, :) = -EEG.data(60, :);
                end
                [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
            end
        end

        % append data from both runs
        if length(ALLEEG) == 2
            EEG = pop_mergeset( ALLEEG, [1  2], 0);
            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'gui','off'); 
        end

        % make epoch with all R-peaks
        EEG = pop_epoch( EEG, {  'Rpeak'  }, [-0.15   1.5], 'newname', 'Rpeak epochs all', 'epochinfo', 'yes');
        % EEG = pop_rmbase( EEG, [-150 -50] ,[]);
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
        
        % calculate heart rate for each epoch
        IBI_tmp = []; RT_interval= [];
        for ep = 1:length(EEG.epoch)
            Rpeaks = find(strcmp(EEG.epoch(ep).eventtype, 'Rpeak'));
            Tpeaks = find(strcmp(EEG.epoch(ep).eventtype, 'Tpeak'));
            if length(Rpeaks) > 1 && length(Tpeaks) >= 1 && Tpeaks(1) < Rpeaks(2)
                IBI_tmp(ep) = EEG.epoch(ep).eventlatency{Rpeaks(2)};
                RT_interval(ep) = EEG.epoch(ep).eventlatency{Tpeaks(1)};
            else
                IBI_tmp(ep) = NaN;
                RT_interval(ep) = NaN;
            end
        end
                
        % find outliers - ectopic beats
        zscore_IBI = (IBI_tmp - mean(IBI_tmp, 'omitnan'))/std(IBI_tmp, 'omitnan');
        IBI_tmp(abs(zscore_IBI) > 4) = NaN;
        zscore_RT_interval = (RT_interval - mean(RT_interval, 'omitnan'))/std(RT_interval, 'omitnan');
        RT_interval(abs(zscore_RT_interval) > 4) = NaN;
        
        
        % Exclude outliers that might refer to ectopic beats or
        % trials where the Rpeak/Tpeak was wrongly detected
        trials2keep = intersect(find(~isnan(IBI_tmp)), find(~isnan(RT_interval)));
        
        EEG = pop_select( EEG, 'trial', trials2keep);
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 4,'gui','off'); 
        
        % exclude outliers from IBI variable
        IBI{t} = IBI_tmp(trials2keep)';
        
        % in order to ensure that each epoch is associated with one
        % T-locked HEP and one IBI we need to remove teh second T event
        % from each R-locked epoch
        Tevent = []; 
        for evnt = 1:length(EEG.event)
            if strcmp(EEG.event(evnt).type,  'Tpeak')
                Tevent = [Tevent; evnt EEG.event(evnt).epoch];
            end
        end
        evnt2del =[];
        for evnt = 2:size(Tevent, 1)
             if Tevent(evnt, 2) == Tevent(evnt-1, 2)
                    evnt2del = [evnt2del; Tevent(evnt, 1)];
             end
        end
        EEG.event(evnt2del) = [];
        
        t_peaks = 0;
        for evn = 1:length(EEG.event)
            if strcmp(EEG.event(evn).type, 'Tpeak')
                t_peaks = t_peaks + 1;
            end
        end
        if t_peaks ~= 0 % in case t peaks were detected
            EEG = pop_epoch( EEG, {  'Tpeak'  }, [-0.2         0.4], 'newname', 'Tpeak epochs b4cue', 'epochinfo', 'yes');
            [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
            % HEP T-locked data
            ECG_T_locked_data = cat(3, ECG_T_locked_data, EEG.data(60, 126:end, :));
            number_trials(t) = size(EEG.data, 3);
        else % in case t peaks were not detected
            % HEP T-locked data
            ECG_T_locked_data = cat(3, ECG_T_locked_data, []);
            number_trials(t) = 0;
        end
    end
    
    if sum(number_trials) > 0    
        % create variable with data electrodes X time frames X trials -
        % starting at 
        % categorical variable with different task conditions
        categorical_variable=[ones(number_trials(1), 1); ones(number_trials(2), 1)*2];

        % create continuous variable with IBI data for each task condition in
        % different columns
        % z-score IBI data
        for t = 1:length(tasks)
            zscore_ibi{t} = zscore(IBI{t});
        end

        continuous_variable = zeros(length(categorical_variable), 2);
        continuous_variable(1:number_trials(1), 1) = zscore_ibi{1};
        continuous_variable(number_trials(1)+1:number_trials(2)+number_trials(1), 2) = zscore_ibi{2};


        LIMO_directory=strcat(pwd, '\AB', num2str(subj_number));
        mkdir(LIMO_directory);

        save(strcat(LIMO_directory, filesep, 'ECG_T_locked_data.mat'), 'ECG_T_locked_data');
        save(strcat(LIMO_directory, filesep, 'categorical_variable.mat'), 'categorical_variable');
        save(strcat(LIMO_directory, filesep, 'continuous_variable.mat'), 'continuous_variable');

        % create LIMO variable
        LIMO.Level                    = 1;

        LIMO.Type = 'Channels';
        LIMO.Analysis = 'Time';

        LIMO.dir                      = LIMO_directory;

        LIMO.data.data_dir            = LIMO_directory;
        LIMO.data.data                = 'ECG_T_locked_data.mat';

        LIMO.data.start               = 51;
        LIMO.data.end                 = 400;
        LIMO.data.trim1               = 1; 
        LIMO.data.trim2               = size(ECG_T_locked_data, 2);

        LIMO.data.timevect            = 51:2:400;

        LIMO.data.sampling_rate       = 500;

        LIMO.data.Cat                 = categorical_variable; 
        LIMO.data.Cont                = continuous_variable; 

        

        LIMO.data.neighbouring_matrix = channeighbstructmat(1, 1);
        LIMO.data.chanlocs = expected_chanlocs(1);

        LIMO.design.fullfactorial     = 0;
        LIMO.design.zscore            = 0;
        LIMO.design.method            = 'OLS'; 
        LIMO.design.type_of_analysis  = 'Mass-univariate';
        LIMO.design.bootstrap         = 0; 
        LIMO.design.tfce              = 0;



        [LIMO.design.X, LIMO.design.nb_conditions, LIMO.design.nb_interactions, LIMO.design.nb_continuous] = limo_design_matrix(ECG_T_locked_data, LIMO, 0);

        LIMO.design.status          = 'to do';
        % LIMO.design.name  = 'paired t-test all electrodes';

        save(strcat(LIMO_directory, filesep, 'LIMO.mat'), 'LIMO');

        % run LIMO first level
        cd(LIMO_directory);
        limo_eeg(4)

        cd ..
    end
end
