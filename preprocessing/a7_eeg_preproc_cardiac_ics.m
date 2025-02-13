% EEG analysis for HEP study - import ECG events
clear; close all
save_dir = 'E:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\cardiac_ICs';
%load eeg data and bad channels data
dirData = dir('E:\ProjectAgeingAuditoryTask\BIDS');      %# Get the data for the current directory
dirIndex = [dirData.isdir];
subDirs = {dirData(dirIndex).name};
validIndex = ~ismember(subDirs,{'.','..'});
% %load files with bad channels
% dirData2 = dir('C:\Users\april\Desktop\erasmus - tesi magistrale\ECG_Toolbox_DEI\ecgSimulator3.0\ecgSimulator3.0\DATA\files_with_BadChannels');      %# Get the data for the current directory
% dirIndex2 = [dirData2.isdir];
% subDirs2 = {dirData2(dirIndex2).name};  
% validIndex2 = ~ismember(subDirs2,{'.','..'});
%load R peak and T peak events
dirData_events = dir('E:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\neweeg_events');      %# Get the data for the current directory
dirIndex3 = [dirData_events.isdir];
subDirs_events = {dirData_events(dirIndex3).name};  
validIndex3 = ~ismember(subDirs_events,{'.','..'});
subDirs_events = subDirs_events(validIndex3);

%load ica files
dirData_ica = dir('E:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\files_with_ica');      %# Get the data for the current directory
dirIndex4 = [dirData_ica.isdir];
fileList_ica = {dirData_ica(~dirIndex4).name}';
%load artifacts
dirData_art_ics = dir('E:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\artifactICs');      %# Get the data for the current directory
dirIndex5 = [dirData_art_ics.isdir];
fileList_a = {dirData_art_ics(~dirIndex5).name}';

% list of participants included
sbj_included = subDirs_events;

% dir with artefact segments
dir_artefact_seg = 'E:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\artefact_segments';

% variable to save independent components identified as cardiac afterfacts
number_cardiac_ics = [];
% column 1 = ics that correlate with ECG
% column 2 = mean absolute correlation coefficient from ICs that correlated
% with ECG
% column 3 = no ic correlates with ECG, 1 was shown as highest
% correlation (max(abs(r)))
% column 4 = cardiac IC already identified as artifact


[ALLEEG EEG CURRENTSET ALLCOM] = eeglab; pop_editoptions( 'option_computeica', 1);
for s = 1:length(sbj_included)% 54%37
%       iDir = find(validIndex3) % loop through participants included and with eeg events calculated

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
   % calculated in E:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\segments_with_artefacts.m
   load([dir_artefact_seg, filesep, sbj_included{s}(5:end), '_artefact_segments.mat'])
   
   
   EKG_data = cell(2, 2);
   
   %% load eeg file
   for k = 1:length(fileList_eeg)
       if contains(fileList_eeg(k), 'eeg.set') && ~contains(fileList_eeg(k), 'passive')
           
           if contains(fileList_eeg(k), 'gonogo')
               task = 2;
           elseif contains(fileList_eeg(k), 'simpleRT')
               task = 1;
           end
           
           if contains(fileList_eeg(k), 'run-1')
               run = 1;
           elseif contains(fileList_eeg(k), 'run-2')
               run = 2;
           end
           

            EEG = pop_loadset('filename',fileList_eeg(k), 'filepath', nextDir_eeg{1});
            [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
            EEG = eeg_checkset( EEG );

           % filter eeg data between .5 and 45 Hz
            EEG = pop_eegfiltnew(EEG, 'locutoff',0.5,'hicutoff',45);
            
            % extract data between start and end triggers (9 and 10)
            start_time = []; end_time = []; cue_times = [];
            for evn = 1:length(EEG.event)
                if EEG.event(evn).type == 9
                    start_time = EEG.event(evn).latency;
                elseif EEG.event(evn).type == 10
                    end_time = EEG.event(evn).latency;
                elseif EEG.event(evn).type == 1
                    cue_times = [cue_times; EEG.event(evn).latency];
                end
            end
            
            if ~isempty(start_time) && ~isempty(end_time)
                epoch_size = (end_time-start_time)/EEG.srate;
                EEG = pop_epoch( EEG, {  '9'  }, [0 epoch_size], 'newname', fileList_eeg{k}(1:end-4), 'epochinfo', 'yes');
            elseif ~isempty(start_time) && isempty(end_time)
                epoch_size = (size(EEG.data, 2)-start_time)/EEG.srate;
                EEG = pop_epoch( EEG, {  '9'  }, [0 epoch_size], 'newname', fileList_eeg{k}(1:end-4), 'epochinfo', 'yes');
            else
                epoch_size = (size(EEG.data, 2)-cue_times(1))/EEG.srate;
                EEG = pop_epoch( EEG, {  '1'  }, [0 epoch_size], 'newname', fileList_eeg{k}(1:end-4), 'epochinfo', 'yes');
            end
            
           % save data from EKG channel into matlab variable - channel 64 (62 after the re-reference)
           EKG_chn = [];
           for chn = 1:length(EEG.chanlocs)
               if strcmp(EEG.chanlocs(chn).labels, 'EKG')
                   EKG_chn = chn;
               end
           end
           
            EKG_data{task, run} = EEG.data(EKG_chn, :);
            EKG_chanlocs = EEG.chanlocs(EKG_chn);
            
            % remove VEO, HEO, EKG and pupil channels
            EEG = pop_select( EEG,'nochannel',{'VEO' 'HEO' 'EKG' 'R-Dia-X-(mm)' 'R-Dia-Y-(mm)'});
            
           % re-reference to linked M1 and M2 (earlobes)
            EEG = pop_reref(EEG, [30 40]);

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
            
            [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
            
       end
   end
   
   % merge dataset
   includ_ds = [];
   for ds = 1:length(ALLEEG)
       if contains(ALLEEG(ds).setname, 'simpleRT_run-1')
           includ_ds(1) = ds;
       elseif contains(ALLEEG(ds).setname, 'simpleRT_run-2')
           includ_ds(2) = ds;
       elseif contains(ALLEEG(ds).setname, 'gonogo_run-1')
           includ_ds(3) = ds;
       elseif contains(ALLEEG(ds).setname, 'gonogo_run-2')
           includ_ds(4) = ds;
       end
   end
   
   EEG = pop_mergeset( ALLEEG, includ_ds, 0);
   [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'setname','task_datasets','gui','off'); 

   % import independent components (ICs) into dataset
   EEG = pop_editset(EEG, 'icachansind', 'ALLEEG(1).icachansind', 'icaweights', 'ALLEEG(1).icaweights', 'icasphere', 'ALLEEG(1).icasphere');

    % identify ICs capturing cardiac artefacts
    
    EKG_data_all = [EKG_data{1, 1}, EKG_data{1, 2}, EKG_data{2, 1}, EKG_data{2, 2}];
    
    clear r p
    for i = 1:size(EEG.icaact, 1)
        [r(i),p(i)] = corr(EKG_data_all', EEG.icaact(i, :)');
    end

    p_thresh = .05; r_thresh = .1;
    % find highly correlated components
    sig_corr = find(p < p_thresh);
    corr_comp = find(abs(r) > r_thresh);



    % variable to save independent components identified as cardiac afterfacts
    cardiac_ics = cell(4, 1);
    % cell 1 = ics that correlate with ECG
    % cell 2 = correlation coefficients from ICs that correlated with ECG
    % cell 3 = no ic correlates with ECG, 1 was shown as highest
    % correlation, ic number and (max(abs(r)))
    % cell 4 = cardiac ICs that had already been identified as artifact
    % cell 5 = cardiac ICs removed from data


    % take comps above threshold
    if ~isempty(intersect(sig_corr,corr_comp))
        artefact_comp = intersect(sig_corr,corr_comp);
        cardiac_ics{1} = artefact_comp;
        cardiac_ics{2} = r(artefact_comp);
        cardiac_ics{3} = [];
        
        number_cardiac_ics(s, 1) = length(artefact_comp);
        number_cardiac_ics(s, 2) = mean(abs(r(artefact_comp)));
        number_cardiac_ics(s, 3) = NaN;
        
    else % take highest one
        artefact_comp = find(abs(r) == max(abs(r)));
        cardiac_ics{1} = [];
        cardiac_ics{2} = [];
        cardiac_ics{3} =  [artefact_comp, max(abs(r))];
        
        number_cardiac_ics(s, 1) = NaN;
        number_cardiac_ics(s, 2) = NaN;
        number_cardiac_ics(s, 3) = max(abs(r));
    end

   % load file with artefact components - include cardiac comp
   str = [dirData_art_ics(1).folder,filesep, [sbj_included{s}(5:end), '_artifact_ICs.mat']];
   artifact_ICs = load(str);
   cardiac_ics{4} = []; number_cardiac_ics(s, 4) = 0;
   for index = 1:length(artefact_comp)
       if any(artifact_ICs.artifact_ICs(:) == artefact_comp(index))
        %do nothing
           cardiac_ics{4} = [cardiac_ics{4}; artefact_comp(index)];
           number_cardiac_ics(s, 4) = number_cardiac_ics(s, 4) + 1;
       else
           artifact_ICs.artifact_ICs(end+1) = artefact_comp(index);
       end
   end
       
   if isempty(cardiac_ics{3})
        cardiac_ics{5} = setdiff(cardiac_ics{1}, cardiac_ics{4});
   else
        cardiac_ics{5} = setdiff(cardiac_ics{3}(1), cardiac_ics{4});
   end
   
   save([save_dir, filesep, sbj_included{s}(5:end), '_cardiac_ICs.mat'], 'cardiac_ics')
end

save([save_dir, filesep, 'number_cardiac_ics.mat'], 'number_cardiac_ics')


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