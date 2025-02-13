% calculate HEP_T_locked 2s interval before cue and 2 s interval after cue
clear; close all;
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

older = [11 12 14 20 21 22 23 32 37 38 41 43 47 48 49 52 55 57 58 63 64 65 67 69 7 70 71 75 8 83 86];
young = [13 15 16 25 26 28 31 33 34 36 4 42 44 45 46 50 51 53 54 56 59 6 62 66 68 72 74 76 78 80 82 84 85 9];
% note: AB21, AB50 and AB55 were excluded from T peak analyses because (s = 8 30 35)
% these were not clearly visible in the ECGs

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
% subjects excluded from T-locked HEP analyses due to difficulty
% identifying AB21, AB50 and AB55
subjects = setdiff(subjects, {'AB21', 'AB50', 'AB55'});

tasks = {'passive', 'simpleRT', 'gonogo'};

HEP_T_locked_all = cell(2, 2); HEP_T_locked_aftercue = cell(2, 2); HEP_T_locked_b4cue = cell(2,2);
number_cue_epochs = cell(2,2); percent_epoch_target = cell(2,2); HEP_target_latency = {};
number_epochs_after_cue = []; number_epochs_b4cue = []; 

IBI_aftercue = cell(2,2); IBI_b4cue = cell(2,2);

yng = 0; old = 0; without_t = [];
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
    
    HEP_T_locked_all{group, 1} = [HEP_T_locked_aftercue{group, 1}; subj_number];
    HEP_T_locked_aftercue{group, 1} = [HEP_T_locked_aftercue{group, 1}; subj_number];
    HEP_T_locked_b4cue{group, 1} = [HEP_T_locked_b4cue{group, 1}; subj_number];
    number_cue_epochs{group, 1} = [number_cue_epochs{group, 1}; subj_number];
    percent_epoch_target{group, 1} = [percent_epoch_target{group, 1}; subj_number];
    
    IBI_aftercue{group, 1} = [IBI_aftercue{group, 1}; subj_number];
    IBI_b4cue{group, 1} = [IBI_b4cue{group, 1}; subj_number];
    
    number_epochs_after_cue(s, 1) = subj_number; 
    number_epochs_b4cue(s, 1) = subj_number;  
       
    for t = 1:length(tasks)
        % clear eeglab
        STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];  

        % load data
        filename = {};
        for f = 1:length(folder_content)
            if contains(folder_content(f).name, subjects(s)) && contains(folder_content(f).name, tasks(t))
                filename = folder_content(f).name;
                EEG = pop_loadset('filename', filename, 'filepath', files_folder );
                % EEG = pop_select( EEG, 'nochannel',{'EKG'});
                [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
            end
        end

        % append data from both runs
        if length(ALLEEG) == 2
            EEG = pop_mergeset( ALLEEG, [1  2], 0);
            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'gui','off'); 
        end
        
        % delete R and T peaks associated with abnormal cardiac cycles or
        % cycles where the cycle was not segmented properly
        outliers_evnt = find_cardiac_cycle_outliers(EEG);
        % delete events associated with outlier cycles
        EEG = pop_selectevent( EEG, 'omitevent', outliers_evnt ,'deleteevents','on');

        EEG = pop_epoch( EEG, {  'Tpeak'  }, [-0.3         0.35], 'newname', 'Tpeak epochs all', 'epochinfo', 'yes');
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
        if group == 1
            HEP_T_locked_all{group, 2}(yng, t, :, :) = mean(EEG.data, 3);
        else
            HEP_T_locked_all{group, 2}(old, t, :, :) = mean(EEG.data, 3);
        end

       
        % make epoch before cue
        for f = 1:length(ALLEEG)
            if contains(ALLEEG(f).setname, 'Merged')
                retrieve_dataset = f;
            end
        end
        
        if t == 1 % passive task
            retrieve_dataset = 1;
        end
        
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 4,'retrieve', retrieve_dataset,'study',0); 
        % delete R and T peaks associated with abnormal cardiac cycles or
        % cycles where the cycle was not segmented properly
        outliers_evnt = find_cardiac_cycle_outliers(EEG);
        % delete events associated with outlier cycles
        EEG = pop_selectevent( EEG, 'omitevent', outliers_evnt ,'deleteevents','on');
        
        EEG = pop_epoch( EEG, {  '1'  }, [-1.65   0], 'newname', 'epochs b4cue', 'epochinfo', 'yes');
        
        number_cue_epoch_b4 = size(EEG.data, 3);
        
        %% calculate IBI before cue
        IBI_b4cue_tmp = []; Rpeaks = [];
        for ep = 1:length(EEG.epoch)
            Rpeaks = find(strcmp(EEG.epoch(ep).eventtype, 'Rpeak'));
            if length(Rpeaks) == 2
                IBI_b4cue_tmp(ep) = EEG.epoch(ep).eventlatency{Rpeaks(2)}-EEG.epoch(ep).eventlatency{Rpeaks(1)};
            elseif length(Rpeaks) == 3
                IBI_b4cue_tmp(ep) = mean([EEG.epoch(ep).eventlatency{Rpeaks(2)}-EEG.epoch(ep).eventlatency{Rpeaks(1)}, ...
                    EEG.epoch(ep).eventlatency{Rpeaks(3)}-EEG.epoch(ep).eventlatency{Rpeaks(2)}]);
            elseif length(Rpeaks) <= 1
                IBI_b4cue_tmp(ep) = NaN;
            end
        end
                
        %% find outliers
        zscore_IBI = (IBI_b4cue_tmp - mean(IBI_b4cue_tmp, 'omitnan'))/std(IBI_b4cue_tmp, 'omitnan');
        IBI_outliers = find(abs(zscore_IBI) > 4);
        IBI_b4cue_tmp(IBI_outliers) = NaN;
        
        EEG = pop_epoch( EEG, {  'Tpeak'  }, [-0.3         0.35], 'newname', 'Tpeak epochs b4cue', 'epochinfo', 'yes');
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
        number_epochs_b4cue(s, t+1) = size(EEG.data, 3);
        if group == 1
            HEP_T_locked_b4cue{group, 2}(yng, t, :, :) = mean(EEG.data, 3);
            IBI_b4cue{group, 2}(yng, t) = mean(IBI_b4cue_tmp, 'omitnan');
        else
            HEP_T_locked_b4cue{group, 2}(old, t, :, :) = mean(EEG.data, 3);
            IBI_b4cue{group, 2}(old, t) = mean(IBI_b4cue_tmp, 'omitnan');
        end

        % make epoch after cue
        for f = 1:length(ALLEEG)
            if contains(ALLEEG(f).setname, 'Merged')
                retrieve_dataset = f;
            end
        end
        
        if t == 1 % passive task
            retrieve_dataset = 1;
        end
        
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 4,'retrieve', retrieve_dataset,'study',0); 
        % delete R and T peaks associated with abnormal cardiac cycles or
        % cycles where the cycle was not segmented properly
        outliers_evnt = find_cardiac_cycle_outliers(EEG);
        % delete events associated with outlier cycles
        EEG = pop_selectevent( EEG, 'omitevent', outliers_evnt ,'deleteevents','on');
        
        EEG = pop_epoch( EEG, {  '1'  }, [0     2], 'newname', 'epochs aftercue', 'epochinfo', 'yes');
        EEG = pop_select( EEG, 'time',[.350 2]);
        
        number_cue_epoch_after = size(EEG.data, 3);
        
        if group == 1
            number_cue_epochs{group, 2}(yng, t, :, :) = [number_cue_epoch_b4, number_cue_epoch_after];
        else
            number_cue_epochs{group, 2}(old, t, :, :) = [number_cue_epoch_b4, number_cue_epoch_after];
        end
        
        % calculate IBI after the cue
        IBI_aftercue_tmp = []; Rpeaks = [];
        for ep = 1:length(EEG.epoch)
            Rpeaks = find(strcmp(EEG.epoch(ep).eventtype, 'Rpeak'));
            if length(Rpeaks) == 2
                IBI_aftercue_tmp(ep) = EEG.epoch(ep).eventlatency{Rpeaks(2)}-EEG.epoch(ep).eventlatency{Rpeaks(1)};
            elseif length(Rpeaks) == 3
                IBI_aftercue_tmp(ep) = mean([EEG.epoch(ep).eventlatency{Rpeaks(2)}-EEG.epoch(ep).eventlatency{Rpeaks(1)}, ...
                    EEG.epoch(ep).eventlatency{Rpeaks(3)}-EEG.epoch(ep).eventlatency{Rpeaks(2)}]);
            elseif isempty(Rpeaks) || length(Rpeaks) == 1
                IBI_aftercue_tmp(ep) = NaN;
            end
        end
                
        % find outliers
        zscore_IBI = (IBI_aftercue_tmp - mean(IBI_aftercue_tmp, 'omitnan'))/std(IBI_aftercue_tmp, 'omitnan');
        IBI_outliers = find(abs(zscore_IBI) > 4);
        IBI_aftercue_tmp(IBI_outliers) = NaN;
        
        
        EEG = pop_epoch( EEG, {  'Tpeak'  }, [-0.3  0.350], 'newname', 'Tpeak epochs aftercue', 'epochinfo', 'yes');
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );

        
        % calculate how many epochs contain the target and the timing of
        % the target and exclude these from the analysis
        if t > 1
            epoch_target = []; %epoch_target_latency = [];
            for epc = 1:length(EEG.epoch)
                if ismember('2', EEG.epoch(epc).eventtype) || ismember('4', EEG.epoch(epc).eventtype)
                    epoch_target = [epoch_target; epc];
%                     epoch_target_latency = [epoch_target_latency; cell2mat(EEG.epoch(epc).eventlatency(ismember(EEG.epoch(epc).eventtype, '2')))];
                end
            end
            
            EEG = pop_select( EEG, 'notrial', epoch_target);
            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 7,'gui','off'); 
            
        end
        
        number_epochs_after_cue(s, t+1) = size(EEG.data, 3);
        if group == 1
            HEP_T_locked_aftercue{group, 2}(yng, t, :, :) = mean(EEG.data, 3);
            IBI_aftercue{group, 2}(yng, t) = mean(IBI_aftercue_tmp, 'omitnan');
        else
            HEP_T_locked_aftercue{group, 2}(old, t, :, :) = mean(EEG.data, 3);
            IBI_aftercue{group, 2}(old, t) = mean(IBI_aftercue_tmp, 'omitnan');
        end

    end       
end

save_dir = 'E:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\TLockedNoBaseline\NoBaseline05HzFilt_NoTarget_cardiacICs';
cd(save_dir)

save number_cue_epochs number_cue_epochs

save HEP_T_locked_aftercue HEP_T_locked_aftercue
save HEP_T_locked_b4cue HEP_T_locked_b4cue
save HEP_T_locked_all HEP_T_locked_all

save number_epochs_after_cue number_epochs_after_cue
save number_epochs_b4cue number_epochs_b4cue

save IBI_aftercue IBI_aftercue
save IBI_b4cue IBI_b4cue

%% plot figures
save_dir = 'E:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\TLockedNoBaseline\NoBaseline05HzFilt_NoTarget_cardiacICs';
cd(save_dir)
chan_dir = 'E:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study';
channels = load([chan_dir filesep 'chanlocs.mat']);
tasks = {'passive', 'simpleRT', 'gonogo'};
load HEP_T_locked_aftercue; load HEP_T_locked_b4cue; load HEP_T_locked_all;

for t = 1:length(tasks)
% 
% figure;
% for p = 1:size(HEP_T_locked_b4cue{2, 2}, 1)
%     plot(squeeze(HEP_T_locked_b4cue{2, 2}(p, t, 7, :)))
%     waitforbuttonpress
% end
% close all

% HEP_T_locked before cue
young_b4cue_mean = squeeze(mean(HEP_T_locked_b4cue{1, 2}(:, t, :, :), 1, 'omitnan'));
older_b4cue_mean = squeeze(mean(HEP_T_locked_b4cue{2, 2}(:, t, :, :), 1, 'omitnan'));

young_b4cue_se = squeeze(std(HEP_T_locked_b4cue{1, 2}(:, t, :, :), [], 1, 'omitnan'))/sqrt(size(HEP_T_locked_b4cue{1, 2}, 1));
older_b4cue_se = squeeze(std(HEP_T_locked_b4cue{2, 2}(:, t, :, :), [], 1, 'omitnan'))/sqrt(size(HEP_T_locked_b4cue{1, 2}, 1));

% figure;
% topoplot(young_b4cue_mean(:, 225), channels.chanlocs); 
% Axis.FontSize = 16; caxis([-2 2]); 
% colorbar;
% % colorbar('Ticks',[0, 2, 4, 6], 'FontSize', 30, 'FontWeight','normal');
% colormap(crameri('hawaii'));
% title(['young b4cue ', tasks{t}] , 'FontSize', 32, 'FontWeight','normal')
% 
% figure;
% topoplot(older_b4cue_mean(:, 225), channels.chanlocs); 
% Axis.FontSize = 16; caxis([-2 2]); 
% colorbar;
% % colorbar('Ticks',[0, 2, 4, 6], 'FontSize', 30, 'FontWeight','normal');
% colormap(crameri('hawaii'));
% title(['older b4cue ', tasks{t}] , 'FontSize', 32, 'FontWeight','normal')


% HEP_T_locked after cue
young_aftercue_mean = squeeze(mean(HEP_T_locked_aftercue{1, 2}(:, t, :, :), 1, 'omitnan'));
older_aftercue_mean = squeeze(mean(HEP_T_locked_aftercue{2, 2}(:, t, :, :), 1, 'omitnan'));

young_aftercue_se = squeeze(std(HEP_T_locked_aftercue{1, 2}(:, t, :, :), [], 1, 'omitnan'))/sqrt(size(HEP_T_locked_aftercue{1, 2}, 1));
older_aftercue_se = squeeze(std(HEP_T_locked_aftercue{2, 2}(:, t, :, :), [], 1, 'omitnan'))/sqrt(size(HEP_T_locked_aftercue{1, 2}, 1));

% figure;
% topoplot(young_aftercue_mean(:, 225), channels.chanlocs); 
% Axis.FontSize = 16; caxis([-2 2]); 
% colorbar;
% % colorbar('Ticks',[0, 2, 4, 6], 'FontSize', 30, 'FontWeight','normal');
% colormap(crameri('hawaii'));
% title(['young aftercue ', tasks{t}] , 'FontSize', 32, 'FontWeight','normal')
% 
% figure;
% topoplot(older_aftercue_mean(:, 225), channels.chanlocs); 
% Axis.FontSize = 16; caxis([-2 2]); 
% colorbar;
% % colorbar('Ticks',[0, 2, 4, 6], 'FontSize', 30, 'FontWeight','normal');
% colormap(crameri('hawaii'));
% title(['older aftercue ', tasks{t}] , 'FontSize', 32, 'FontWeight','normal')

% HEP_T_locked all
young_all_mean = squeeze(mean(HEP_T_locked_all{1, 2}(:, t, :, :), 1, 'omitnan'));
older_all_mean = squeeze(mean(HEP_T_locked_all{2, 2}(:, t, :, :), 1, 'omitnan'));

young_all_se = squeeze(std(HEP_T_locked_all{1, 2}(:, t, :, :), [], 1, 'omitnan'))/sqrt(size(HEP_T_locked_all{1, 2}, 1));
older_all_se = squeeze(std(HEP_T_locked_all{2, 2}(:, t, :, :), [], 1, 'omitnan'))/sqrt(size(HEP_T_locked_all{1, 2}, 1));


% figure;
% topoplot(young_all_mean(:, 225), channels.chanlocs); 
% Axis.FontSize = 16; caxis([-2 2]); 
% colorbar;
% % colorbar('Ticks',[0, 2, 4, 6], 'FontSize', 30, 'FontWeight','normal');
% colormap(crameri('hawaii'));
% title(['young ', tasks{t}] , 'FontSize', 32, 'FontWeight','normal')
% 
% figure;
% topoplot(older_all_mean(:, 225), channels.chanlocs); 
% Axis.FontSize = 16; caxis([-2 2]); 
% colorbar;
% % colorbar('Ticks',[0, 2, 4, 6], 'FontSize', 30, 'FontWeight','normal');
% colormap(crameri('hawaii'));
% title(['older ', tasks{t}] , 'FontSize', 32, 'FontWeight','normal')


% channels = F3,Fz,F4,C3,Cz,C4,P3,Pz,P4
channels_num = 34;% [7, 16, 25, 34];

    for p = channels_num

        x_axis=-0.299:0.002:0.349; % in seconds
        figure;
        plot(x_axis, young_b4cue_mean(p, :), 'color', [0 0.4470 0.7410],  'LineWidth',1.5) %  (light blue)
        hold on
        jbfill(x_axis, [young_b4cue_mean(p, :) + young_b4cue_se(p, :)], [young_b4cue_mean(p, :) - young_b4cue_se(p, :)], [0 0.4470 0.7410],[0 0.4470 0.7410], 1, .2)
        hold on
        plot(x_axis, young_aftercue_mean(p, :), 'color', 'k',  'LineWidth',1.5)
        hold on
        jbfill(x_axis, [young_aftercue_mean(p, :) + young_aftercue_se(p, :)], [young_aftercue_mean(p, :) - young_aftercue_se(p, :)], 'k', 'k', 1, .2)
        hold on
        plot(x_axis, zeros(length(x_axis), 1))
        hold off
        box off
        ax = gca;
%         if t == 1
%             axis([-inf inf -2 2])
%         else
%             axis([-inf inf -3 1])
%         end
        ax.LineWidth = 2.5;
        ax.FontSize = 28;
        ax.FontName = 'Arial';
        xlabel('Time (s)', 'FontSize', 32, 'FontWeight','normal')
        ylabel('Amplitude (\muV)', 'FontSize', 32, 'FontWeight','normal')
        name_channels =  channels.chanlocs(p).labels;
        title(['young ', tasks{t}, name_channels] , 'FontSize', 32, 'FontWeight','normal')
%         basefilename2 = sprintf('Tepochplot_gng-sRT_young_cut_%d.fig', p);
%         savefig(basefilename2);


        figure;
        plot(x_axis, older_b4cue_mean(p, :), 'color', [0 0.4470 0.7410],  'LineWidth',1.5) %plot mean gng (light blue)
        hold on
        jbfill(x_axis, [older_b4cue_mean(p, :) + older_b4cue_se(p, :)], [older_b4cue_mean(p, :) - older_b4cue_se(p, :)], [0 0.4470 0.7410],[0 0.4470 0.7410], 1, .2)
        hold on
        plot(x_axis, older_aftercue_mean(p, :), 'color', 'k',  'LineWidth',1.5) %plot mean gng (light blue)
        hold on
        jbfill(x_axis, [older_aftercue_mean(p, :) + older_aftercue_se(p, :)], [older_aftercue_mean(p, :) - older_aftercue_se(p, :)], 'k', 'k', 1, .2)
        hold on
        plot(x_axis, zeros(length(x_axis), 1))
        hold off
        box off
        ax = gca;
        % if t == 1
        %     axis([-inf inf -2 2])
        % else
        %     axis([-inf inf -3 1])
        % end
        ax.LineWidth = 2.5;
        ax.FontSize = 28;
        ax.FontName = 'Arial';
        xlabel('Time (s)', 'FontSize', 32, 'FontWeight','normal')
        ylabel('Amplitude (\muV)', 'FontSize', 32, 'FontWeight','normal')
        name_channels =  channels.chanlocs(p).labels;
        title(['older ', tasks{t}, name_channels] , 'FontSize', 32, 'FontWeight','normal')
        
        
        figure;
        plot(x_axis, young_all_mean(p, :), 'color', 'k',  'LineWidth',1.5) %plot(light blue)
        hold on
        jbfill(x_axis, [young_all_mean(p, :) + young_all_se(p, :)], [young_all_mean(p, :) - young_all_se(p, :)], 'k', 'k', 1, .2)
        hold on
        plot(x_axis, older_all_mean(p, :), '--r',  'LineWidth',1.5) 
        hold on
        jbfill(x_axis, [older_all_mean(p, :) + older_all_se(p, :)], [older_all_mean(p, :) - older_all_se(p, :)], 'r', 'r', 1, .2)
        hold on
        plot(x_axis, zeros(length(x_axis), 1))
        hold off
        box off
        ax = gca;
%         if t == 1
%             axis([-inf inf -2 2])
%         else
%             axis([-inf inf -3 1])
%         end
        ax.LineWidth = 2.5;
        ax.FontSize = 28;
        ax.FontName = 'Arial';
        xlabel('Time (s)', 'FontSize', 32, 'FontWeight','normal')
        ylabel('Amplitude (\muV)', 'FontSize', 32, 'FontWeight','normal')
        name_channels =  channels.chanlocs(p).labels;
        title(['HEP_T_locked ', tasks{t}, name_channels] , 'FontSize', 32, 'FontWeight','normal')

    end

end

%% mean and sd of number of epochs used - before the cue
epochs_young = []; epochs_older = [];
for s = 1:size(number_epochs_b4cue, 1)
    if ismember(number_epochs_b4cue(s, 1), young)
        epochs_young = [epochs_young; number_epochs_b4cue(s, 3:4)];
    else
        epochs_older = [epochs_older; number_epochs_b4cue(s, 3:4)];
    end
end

simpleRT_mean1 = mean(epochs_young(:, 1), 1)
simpleRT_mean2 = mean(epochs_older(:, 1), 1)

simpleRT_sd1 = std(epochs_young(:, 1), [], 1)
simpleRT_sd2 = std(epochs_older(:, 1), [], 1)

gng_mean1 = mean(epochs_young(:, 2), 1)
gng_mean2 = mean(epochs_older(:, 2), 1)

gng_sd1 = std(epochs_young(:, 2), [], 1)
gng_sd2 = std(epochs_older(:, 2), [], 1)


%% epochs after the cue
clc
epochs_young = []; epochs_older = [];
for s = 1:size(number_epochs_after_cue, 1)
    if ismember(number_epochs_after_cue(s, 1), young)
        epochs_young = [epochs_young; number_epochs_after_cue(s, 3:4)];
    else
        epochs_older = [epochs_older; number_epochs_after_cue(s, 3:4)];
    end
end

simpleRT_mean1 = mean(epochs_young(:, 1), 1)
simpleRT_mean2 = mean(epochs_older(:, 1), 1)

simpleRT_sd1 = std(epochs_young(:, 1), [], 1)
simpleRT_sd2 = std(epochs_older(:, 1), [], 1)

gng_mean1 = mean(epochs_young(:, 2), 1)
gng_mean2 = mean(epochs_older(:, 2), 1)

gng_sd1 = std(epochs_young(:, 2), [], 1)
gng_sd2 = std(epochs_older(:, 2), [], 1)


%% calculate number of epochs before and after the cue

% number_cue_epochs{group, 2}(yng, t, :, :) = [number_cue_epoch_b4, number_cue_epoch_after];

% before after cue

simpleRT_mean1 = squeeze(mean(number_cue_epochs{1, 2}(:, 2, :, :), 1))
simpleRT_mean2 = squeeze(mean(number_cue_epochs{2, 2}(:, 2, :, :), 1))

simpleRT_sd1 = std(number_cue_epochs{1, 2}(:, 2, :, :), [], 1)
simpleRT_sd2 = std(number_cue_epochs{2, 2}(:, 2, :, :), [], 1)

gng_mean1 = squeeze(mean(number_cue_epochs{1, 2}(:, 3, :, :), 1))
gng_mean2 = squeeze(mean(number_cue_epochs{2, 2}(:, 3, :, :), 1))

gng_sd1 = std(number_cue_epochs{1, 2}(:, 3, :, :), [], 1)
gng_sd2 = std(number_cue_epochs{2, 2}(:, 3, :, :), [], 1)

%% calculate IBI after and before cue
save_dir = 'E:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\TLockedNoBaseline\NoBaseline05HzFilt_NoTarget_cardiacICs';
cd(save_dir)
% save data in excel file to perform stats in SPSS
load IBI_aftercue; load IBI_b4cue;
% create table
% participant, group, IBI_b4_passive,IBI_after_passive, IBI_b4_simpleRT,IBI_after_simpleRT, IBI_b4_gng,IBI_after_gng, 
% IBI_avg_std{group, 2}(old, t, :) = [mean(IBI), std(IBI)];
IBI_table = table([IBI_b4cue{1, 1};IBI_b4cue{2, 1}],...
    [ones(length(IBI_b4cue{1, 1}), 1); ones(length(IBI_b4cue{2, 1}), 1)*2], ...
    [IBI_b4cue{1, 2}(:, 1); IBI_b4cue{2, 2}(:, 1)], ...
    [IBI_aftercue{1, 2}(:, 1); IBI_aftercue{2, 2}(:, 1)], ...
    [IBI_b4cue{1, 2}(:, 2); IBI_b4cue{2, 2}(:, 2)], ...
    [IBI_aftercue{1, 2}(:, 2); IBI_aftercue{2, 2}(:, 2)], ...
    [IBI_b4cue{1, 2}(:, 3); IBI_b4cue{2, 2}(:, 3)], ...
    [IBI_aftercue{1, 2}(:, 3); IBI_aftercue{2, 2}(:, 3)], ...
    'VariableNames',{'participant', 'group', ' IBI_b4_passive' , 'IBI_after_passive',...
    'IBI_b4_simpleRT', 'IBI_after_simpleRT', 'IBI_b4_gng', 'IBI_after_gng'});

filename = 'IBI_b4_after_cue.xlsx';
writetable(IBI_table,filename)




mean_young_b4cue = mean(IBI_b4cue)
