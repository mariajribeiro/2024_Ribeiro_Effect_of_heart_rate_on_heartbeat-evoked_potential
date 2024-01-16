% find data segments that were manually rejected
% EEG.event = boundary
% duration in data points
% latency - start of segment in data points
clear; close all
sbj_included = 43;%[ 4,6,7,8,9,11,12,13,14,15,16,22,23,26,28,31,32,33,34,36,41,46,47,48,49 ...
% 51,52,53,54,55,57,58,62,63,64,65,66,68,70,72,75,78,80,83,85,86 ...
% 20,21,25,37,38,42,44,45,50,56,59,67,69,71,74,76,82,84];

save_dir = 'D:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\artefact_segments';
data_folder = 'G:\ProjectAgingNeuromodulation\AuditoryResearch\AuditoryTask_EyeTracking_EEG_lab94';
subjct_folder = dir(data_folder);
task={'W1', 'D1', 'D2', 'G1', 'G2'};

[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

for s = 1:length(sbj_included)
    
    if ismember(sbj_included(s), [14, 15])
        task={'W2', 'D1', 'D2', 'G1', 'G2'};
    else
        task={'W1', 'D1', 'D2', 'G1', 'G2'};
    end
    
    filepath = [data_folder, filesep, ['AB', num2str(sbj_included(s))], filesep, 'EEG'];
    artefact_segments = {};
    for t = 1:length(task)
        % clear eeglab
        STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
        
        filename = ['AB', num2str(sbj_included(s)), '_', task{t}, '_Filt0_1_35Hz_Pupil_ManualArtRej.set'];
        EEG = pop_loadset(filename, filepath);        % load dataset with artifact segments removed

        boundary_start = []; boundary_duration = [];
        for evn = 1:length(EEG.event)
            if strcmp(EEG.event(evn).type, 'boundary')
                boundary_start = [boundary_start; round(EEG.event(evn).latency)];
                boundary_duration = [boundary_duration; round(EEG.event(evn).duration)];
            end
        end

        boundary_start_end = [boundary_start, boundary_start+boundary_duration-1];
        
        artefact_segments{t} = boundary_start_end;
        
%         EEG = eeg_eegrej( EEG, [671 1099;1924 2307]);
    end
    
    save([save_dir, filesep, 'AB', num2str(sbj_included(s)), '_artefact_segments'], 'artefact_segments');
end



% for seg = 1:size(artefact_segments{1}, 1)
%     EEG = eeg_eegrej( EEG, artefact_segments{1}(seg, :));
% end