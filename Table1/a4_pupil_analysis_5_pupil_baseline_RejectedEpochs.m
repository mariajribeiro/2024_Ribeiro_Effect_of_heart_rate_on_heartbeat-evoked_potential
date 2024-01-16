% determining baseline pupil diameter - epochs rejected
clear
% subjects with pupil data
young_pupil=[4	9	10	13	15	16	25	26	28	31	33	34	36	42	44	45	46 50 51 53 54 56	59	62 66	68	72	74	76  78  80 81   82  84  85];
older_pupil=[7  8	11	12	14	17	19	20	21	22	23	32	35	37	41	43	47 48 49 52 55	57	58	60	61	63	64  65	67	69	70	71	73	75  77 79   83  86];
participants=[young_pupil older_pupil];
task={'W1', 'D1', 'D2', 'G1', 'G2'};
% open eeglab
[ALLEEG , ~, CURRENTSET ALLCOM] = eeglab;

for p=41%participants;
    clearvars -except p participants task
    
directory=strcat('M:\ProjectAgingNeuromodulation\AuditoryResearch\AuditoryTask_EyeTracking_EEG_lab94\AB', num2str(p), '\ET\');
cd(directory);

load rejected_epochs

% make variable with epochs rejected including artifact epochs and errors and trials after errors per task run
    Pupil_TotalRejectEpochs{1}=rejected_epochs{1}';
            
    for t=2:5;
                    % load behavior and calculate error trials and trials after error 
                    behavioural_dir=strcat('M:\ProjectAgingNeuromodulation\AuditoryResearch\AuditoryTask_EyeTracking_EEG_lab94\AB', num2str(p), '\ET\',  task{t},'_behavior\');
                    % error_trials - all errors, response2cue, multiple responses, misses, slow responses
                    load([behavioural_dir, 'misses.mat']);
                    load([behavioural_dir, 'multiple_responses.mat']);
                    load([behavioural_dir, 'response2cue.mat']);
                    load([behavioural_dir, 'slow_responses.mat']);

            if    t==2 || t==3;
                error_trials=cat(1, misses, multiple_responses, response2cue, slow_responses);
                Pupil_TotalRejectEpochs{t}=sortrows(cat(1, rejected_epochs{t}', [error_trials; error_trials+1]));
            elseif t==4 || t==5;
                load([behavioural_dir, 'response2nogo.mat']);
                error_trials=cat(1, misses, multiple_responses, response2cue, slow_responses, response2nogo);  
                Pupil_TotalRejectEpochs{t}=sortrows(cat(1, rejected_epochs{t}', [error_trials; error_trials+1]));
            end
    end
    
    % in the case of the last trial being an error need to delete the error
    % trial +1 as it would be outside the existign trials
    for t=2:5;
      if isempty(setdiff(61, Pupil_TotalRejectEpochs{t}));
        Pupil_TotalRejectEpochs{t}(Pupil_TotalRejectEpochs{t}==61)=[];
      end
    end
    
    save Pupil_TotalRejectEpochs Pupil_TotalRejectEpochs;
    
    % load raw pupil data per trials and save variable with pupil baseline
    % data
    
    for t=1:5
        if t==1
            filename=strcat('AB', num2str(p), 'W_CueLocked_RawPupilDiameter');
            load(filename);
            PupilBaselinePerTaskRun{t}=mean(W_CueLocked_RawPupilDiameter(192:240, :), 1);
        else
            %clear eeglab
            STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
            filename=strcat('AB', num2str(p), '_', task(t), '_epochs_BadEpochsAndErrorsRemoved.set');
            EEG = pop_loadset('filename',filename,'filepath',directory);
            % save data of cue locked raw pupil diameter in mm
            clear PupilBaseline_Data;
            PupilBaseline_Data=squeeze(EEG.data(3, :,:));
            PupilBaselinePerTaskRun{t}=mean(PupilBaseline_Data(192:240, :), 1);
        end
    end
    save PupilBaselinePerTaskRun PupilBaselinePerTaskRun;
    
end


    