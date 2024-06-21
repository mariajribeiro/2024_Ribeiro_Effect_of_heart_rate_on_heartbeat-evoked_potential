% EEG analysis for HEP study
clear; close all
save_dir = 'M:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\eeg_derivatives_05Hzfilter\random_R_T_events';
%load eeg data and bad channels data
dirData = dir('M:\ProjectAgeingAuditoryTask\BIDS');      %# Get the data for the current directory
dirIndex = [dirData.isdir];
subDirs = {dirData(dirIndex).name};
validIndex = ~ismember(subDirs,{'.','..'});
% %load files with bad channels
% dirData2 = dir('C:\Users\april\Desktop\erasmus - tesi magistrale\ECG_Toolbox_DEI\ecgSimulator3.0\ecgSimulator3.0\DATA\files_with_BadChannels');      %# Get the data for the current directory
% dirIndex2 = [dirData2.isdir];
% subDirs2 = {dirData2(dirIndex2).name};  
% validIndex2 = ~ismember(subDirs2,{'.','..'});
%load R peak and T peak events randomized
dirData_events = dir('M:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\neweeg_events_random');
dirIndex3 = [dirData_events.isdir];
subDirs_events = {dirData_events(dirIndex3).name};  
validIndex3 = ~ismember(subDirs_events,{'.','..'});
subDirs_events = subDirs_events(validIndex3);

%load ica files
dirData_ica = dir('M:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\files_with_ica');      %# Get the data for the current directory
dirIndex4 = [dirData_ica.isdir];
fileList_ica = {dirData_ica(~dirIndex4).name}';
%load artifacts
dirData_art_ics = dir('M:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\artifactICs');      %# Get the data for the current directory
dirIndex5 = [dirData_art_ics.isdir];
fileList_a = {dirData_art_ics(~dirIndex5).name}';

% list of participants included
sbj_included = subDirs_events;

% dir with artefact segments
dir_artefact_seg = 'M:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\artefact_segments';

[ALLEEG EEG CURRENTSET ALLCOM] = eeglab; pop_editoptions( 'option_computeica', 1);
for s = 1:length(sbj_included)

    %subfolders of eeg data
    nextDir_eeg = fullfile(dirData(1).folder, sbj_included(s), 'eeg');    
    nextdirData = dir(nextDir_eeg{1});
    nextdirIndex = [nextdirData.isdir];
    fileList_eeg = {nextdirData(~nextdirIndex).name}';
%     %subfolders of files with bad channels
%     nextDir2 = fullfile('C:\Users\april\Desktop\erasmus - tesi magistrale\ECG_Toolbox_DEI\ecgSimulator3.0\ecgSimulator3.0\DATA\files_with_BadChannels',subDirs2{iDir});    
%     nextdirData2 = dir(nextDir2);
%     nextdirIndex2 = [nextdirData2.isdir];
%     fileList2 = {nextdirData2(~nextdirIndex2).name}';
    %subfolders of events
    nextDir3 = fullfile(dirData_events(1).folder, sbj_included(s));    
    nextdirData3 = dir(char(nextDir3));
    nextdirIndex3 = [nextdirData3.isdir];
    fileList_events = {nextdirData3(~nextdirIndex3).name}';
    % clear eeglab
    STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];   
    % load file with ICA weights
    for ic = 1:length(fileList_ica)
        if contains(fileList_ica(ic), sbj_included{s}(5:end)) && contains(fileList_ica(ic), '.set') 
            EEG = pop_loadset('filename', fileList_ica(ic),'filepath',dirData_ica(1).folder);
        end
    end
   [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG);
   
   % load variable with artefact segments to remove later
   % calculated in M:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\segments_with_artefacts.m
   load([dir_artefact_seg, filesep, sbj_included{s}(5:end), '_artefact_segments.mat'])
   
   %% load eeg file
   for k = 1:length(fileList_eeg)
       if contains(fileList_eeg(k), 'eeg.set')
            EEG = pop_loadset('filename',fileList_eeg(k), 'filepath', nextDir_eeg{1});
            [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
            EEG = eeg_checkset( EEG );

           % filter eeg data between .5 and 45 Hz
            EEG = pop_eegfiltnew(EEG, 'locutoff',0.5,'hicutoff',45);

           % re-reference to linked M1 and M2 (earlobes)
            EEG = pop_reref(EEG, [30 40]);

           % save data from EKG channel into matlab variable - channel 64 (62 after the re-reference)
            EKG_data = EEG.data(62, :);
            EKG_chanlocs = EEG.chanlocs(62);

           % remove VEO, HEO, EKG and pupil channels
            EEG = pop_select( EEG,'nochannel',{'VEO' 'HEO' 'EKG' 'R-Dia-X-(mm)' 'R-Dia-Y-(mm)'});
            
            % create variable with chanlocs with all channels to interpolate bad channels later
            chanlocs = EEG.chanlocs;

           % remove bad channels if necessary
           file_eeg = char(fileList_eeg(k));
           start_str = strfind(file_eeg , 'task');
           end_str = strfind(file_eeg , 'eeg');

           for fl = 1:length(fileList_eeg)
               if contains(fileList_eeg{fl}, file_eeg(start_str:end_str-2)) && contains(char(fileList_eeg(fl)), 'eeg.json')
                    fullfilename_badchannels = fullfile(nextDir_eeg, fileList_eeg(fl));
               end
           end
            clear BadChannels
            fid = fopen(fullfilename_badchannels{1});
            raw = fread(fid,inf);
            str = char(raw');
            fclose(fid);
            val = jsondecode(str);
            if isfield(val,'SubjectArtefactDescription')
                BadChannels = val.SubjectArtefactDescription;
                BadChannels = BadChannels(16:end-1);
            end
            B = exist('BadChannels');
            num_BadChannel = [];
            if B == 1
                singleBadChannel = split(BadChannels, ' ');
                num_BadChannel = zeros(1,length(singleBadChannel));
                for l = 1:length(singleBadChannel)
                   for j = 1:EEG.nbchan
                       num = EEG.chanlocs(j).labels;
                       t = strcmp(num,singleBadChannel(l));
                       if t == 1
                          num_BadChannel(l)=j;
                       end
                   end
                end
                EEG = pop_select( EEG,'nochannel',num_BadChannel);
            end

            % import independent components (ICs) into dataset
            EEG = pop_editset(EEG, 'icachansind', 'ALLEEG(1).icachansind', 'icaweights', 'ALLEEG(1).icaweights', 'icasphere', 'ALLEEG(1).icasphere');

            % identify ICs capturing cardiac artefacts
            clear r p
            for i = 1:size(EEG.icaact, 1)
                [r(i),p(i)] = corr(EKG_data', EEG.icaact(i, :)');
            end

            p_thresh = .05; r_thresh = .1;
            % find highly correlated components
            sig_corr = find(p < p_thresh);
            corr_comp = find(abs(r) > r_thresh);

            % take comps above threshold
            if ~isempty(intersect(sig_corr,corr_comp))
                artefact_comp = intersect(sig_corr,corr_comp);
            else % take highest one
                artefact_comp = find(abs(r) == max(abs(r)));
            end

           % load file with artefact components - include cardiac comp
           str = [dirData_art_ics(1).folder,filesep, [sbj_included{s}(5:end), '_artifact_ICs.mat']];
           artifact_ICs = load(str);
           for index = 1:length(artefact_comp)
               if any(artifact_ICs.artifact_ICs(:) == artefact_comp(index))
                %do nothing
               else
                    artifact_ICs.artifact_ICs(end+1) = artefact_comp(index);
               end
           end
           EEG = pop_subcomp( EEG, artifact_ICs.artifact_ICs, 0);
           
           if B == 1 % if bad channels were removed
            % interpolate bad channels
            [EEG] = interpol(EEG, chanlocs);
           end
           
           % import EKG channel data back into file
            EEG.data(end+1,:) = EKG_data;
            EEG.nbchan = size(EEG.data,1);
            if ~isempty(EEG.chanlocs)
                EEG.chanlocs(end+1) = EKG_chanlocs;
            end

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

           %% create epoch between event 9 and event 10
            N_samples = end_latency-start_latency;
%             end_epoch = N_samples*0.004;
           end_epoch = N_samples/EEG.srate;
           [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
           EEG = eeg_checkset( EEG );
           name=fileList_eeg{k};
           if start_latency == 0
               EEG = pop_epoch( EEG, [0      end_epoch], 'newname', name, 'epochinfo', 'yes');
               [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'gui','off'); 
               EEG = eeg_checkset( EEG );
           else
               EEG = pop_epoch( EEG, {  '9'  }, [0      end_epoch], 'newname', name, 'epochinfo', 'yes');
               [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'gui','off'); 
               EEG = eeg_checkset( EEG );
           end

           %% import events R peaks etc.. into data
           file_eeg = char(fileList_eeg(k));
           start_str = strfind(file_eeg , 'task');
           end_str = strfind(file_eeg , 'eeg');
           for fl = 1:length(fileList_events)
               if contains(char(fileList_events(fl)), [file_eeg(start_str:end_str-4), file_eeg(end_str-2)])
                    fullfilename_events = fullfile(nextDir3,fileList_events(fl));
               end
           end
           
           [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
           EEG = eeg_checkset( EEG );
           EEG = pop_importevent( EEG, 'append','yes','event',fullfilename_events{1},'fields',{'latency','type'},'skipline',1,'timeunit',1);
           
%             EEG = pop_importevent( EEG, 'append','no','event',fullfilename_events{1},'fields',{'latency','type'},'skipline',1,'timeunit',1);
          
           [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
           EEG = eeg_checkset( EEG );
           
           %%  remove segments with artefacts using artefact_segments
           % These files are structures with 5 matrices in:
           % artefact_segments{1} – passive task
           % artefact_segments{2} – simple RT run 1
           % artefact_segments{3} – simple RT run 2
           % artefact_segments{4} – go/no-go run 1
           % artefact_segments{5} – go/no-go run 2
           if contains(file_eeg, 'passive')
               art_seg = artefact_segments{1};
           elseif contains(file_eeg, 'simpleRT_run-1')
               art_seg = artefact_segments{2};
           elseif contains(file_eeg, 'simpleRT_run-2')
                art_seg = artefact_segments{3};
           elseif contains(file_eeg, 'gonogo_run-1')
               art_seg = artefact_segments{4};
           elseif contains(file_eeg, 'gonogo_run-2')
                art_seg = artefact_segments{5};   
           end
           
           sgm2delete = [];  duration2del = []; duration2deltruncated = []; sgm_at_beginning = [];
           for sgm = 1:size(art_seg, 1)
               if art_seg(sgm, 2) <= start_latency
                   sgm2delete = [sgm2delete, sgm];
                   duration2del = [duration2del, art_seg(sgm, 2)-art_seg(sgm, 1)];
               elseif art_seg(sgm, 1) < start_latency && art_seg(sgm, 2) > start_latency
                   sgm_at_beginning = sgm;
                   duration2deltruncated = [duration2deltruncated, start_latency - art_seg(sgm, 1)];
%                    art_seg(sgm, 1) = start_latency;
               end
           end
           
           art_seg = art_seg - (start_latency - sum(duration2del) - sum(duration2deltruncated)) + 1; % to compensate for starting the recording at event 9
           if ~isempty(sgm_at_beginning)
            art_seg(sgm_at_beginning, 1) = 1;
           end
           art_seg(sgm2delete, :) = [];
           
           for seg = 1:size(art_seg, 1)
                EEG = eeg_eegrej( EEG, art_seg(seg, :));
           end
 %%               
%            %remove EKG
%            EEG = pop_select( EEG, 'nochannel',{'EKG'});
%            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 3,'gui','off'); 
           EEG = eeg_checkset( EEG );
           basefilename = sprintf('eeg_random_Rpeak_Tpeak_%s', fileList_eeg{k});
           EEG = pop_saveset( EEG, 'filename',basefilename,'filepath',save_dir);
           [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
       end  
   end
end


%% function to interpolate removed channels - needs chanloc variable with all channels
function [EEG] = interpol( EEG, chanlocs )
    %% interpolation
    if nargin < 2
        load('chanlocs.mat')
    end

    chans_eeg = [];
    for i=1:length(EEG.chanlocs)
        chans_eeg = [ chans_eeg {EEG.chanlocs(i).labels} ];
    end

    idxs = [];
    for i=1:length(chanlocs)
        index = find(ismember(chans_eeg, chanlocs(i).labels) == 1, 1);
        if isempty(index)
            idxs = [idxs i];
        end
    end

    EEG = pop_interp(EEG, chanlocs(idxs), 'spherical');

    % reorder
    chans_eeg = [];
    for c=1:length(EEG.chanlocs)
        chans_eeg = [ chans_eeg {EEG.chanlocs(c).labels} ];
    end

    idxs = [];
    for c=1:length(chanlocs)
        index = find(ismember(chans_eeg, chanlocs(c).labels) == 1, 1);
        idxs = [idxs index];
    end

   EEG.data = EEG.data(idxs,:,:);
   EEG.chanlocs = EEG.chanlocs(idxs);

   indcomps = [];
   for compidx = 1:length(EEG.icachansind)
       indcomps = [indcomps find(EEG.icachansind(compidx) == idxs)];
   end
   EEG.icachansind = indcomps;

end