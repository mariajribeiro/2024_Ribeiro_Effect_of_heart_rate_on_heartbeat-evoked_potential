% create R and T events for subject AB74 - ECG with unsual shape not
% recognized by toolbox
% latency in seconds and event type Rpeak or Tpeak
% Maria Ribeiro
clear; close all
% Data from directory
eeg_data_dir = 'D:\ProjectAgeingAuditoryTask\BIDS\sub-AB74\eeg';
dirData = dir(eeg_data_dir);      %# Get the data for the current directory

% open eeglab
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;


for f = 1:length(dirData)
    if contains(dirData(f).name, 'eeg.set')

        STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];
        EEG = pop_loadset('filename',dirData(f).name,'filepath',eeg_data_dir);
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );

        clear start_latency end_latency

        % remove data before event 9 if exists
        for evn = 1:length(EEG.event)
            if EEG.event(evn).type == 9
               start_latency = EEG.event(evn).latency;
            elseif EEG.event(evn).type == 10
                   end_latency = EEG.event(evn).latency;
            end
        end
        A = exist ('start_latency');
        if A == 0
           start_latency = 0;
        end

        %create epoch between event 9 and event 10
        N_samples = end_latency-start_latency;
        end_epoch = N_samples/EEG.srate;
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
        EEG = eeg_checkset( EEG );
        if start_latency == 0
           EEG = pop_epoch( EEG, [0      end_epoch], 'newname', 'start_events', 'epochinfo', 'yes');
           [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'gui','off'); 
           EEG = eeg_checkset( EEG );
        else
           EEG = pop_epoch( EEG, {  '9'  }, [0      end_epoch], 'newname', 'start_events', 'epochinfo', 'yes');
           [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'gui','off'); 
           EEG = eeg_checkset( EEG );
        end


        % select EKG and resample 
        EEG = pop_select( EEG, 'channel',{'EKG'});
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'gui','off'); 
        EEG = eeg_checkset( EEG );
        EEG = pop_eegfiltnew(EEG, 'locutoff',0.1,'hicutoff',35,'plotfreqz',1);
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 3,'setname','filtered','gui','off'); 
        
        % find R and T peaks as local maxima and minima
        data_ecg = squeeze(EEG.data);
        [pks, locs] = findpeaks(-data_ecg, 'MinPeakProminence',50, 'MinPeakDistance', 200, 'MaxPeakWidth', 20);
        Rpks = pks;
        Rlocs = locs;
        
        [pks, locs] = findpeaks(-data_ecg, 'MinPeakProminence',50, 'MinPeakDistance', 200, 'MinPeakWidth', 30);
        
        Tpks = pks;
        Tlocs = locs;
        
        figure; plot(-data_ecg)
        hold on
        plot(Tlocs, Tpks, 'o')
        
        Tlatency = Tlocs/EEG.srate;
        Rlatency = Rlocs/EEG.srate;
        
%         RR = diff(Rlatency);
%         zscore_RR = zscore(RR)
%         outliers = find(zscore_RR < -2);
%         
%         Rlatency(outliers) = [];
        
        % save R peak events as txt file
        filename = ['newRpeak_Tpeak_', dirData(f).name(1:end-10), dirData(f).name(end-8),'.txt'];
        filepath = 'D:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\neweeg_events\sub-AB74';
        fileID = fopen([filepath, filesep, filename],'wt');
        fprintf(fileID,'%s %s\n','latency','type');
        for i = 1:length(Rlatency)
        fprintf(fileID,'%f %s\n', Rlatency(i),'Rpeak');
        end
        for i = 1:length( Tlatency)
        fprintf(fileID,'%f %s\n',  Tlatency(i),'Tpeak');
        end
        fclose(fileID);

    end
end