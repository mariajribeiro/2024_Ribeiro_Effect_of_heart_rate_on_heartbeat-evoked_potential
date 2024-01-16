% calculate HEP_T_locked 2s interval before cue and 2 s interval after cue
clear; close all;
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

older = [11 12 14 20 21 22 23 32 37 38 41 43 47 48 49 52 55 57 58 63 64 65 67 69 7 70 71 75 8 83 86];
young = [13 15 16 25 26 28 31 33 34 36 4 42 44 45 46 50 51 53 54 56 59 6 62 66 68 72 74 76 78 80 82 84 85 9];
% note: AB21, AB50 and AB55 were excluded from T peak analyses because (s = 8 30 35)
% these were not clearly visible in the ECGs

% all files
files_folder = 'D:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\eeg_derivatives';
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
number_epochs_after_cue = []; number_epochs_b4cue = []; 
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
                EEG = pop_select( EEG, 'nochannel',{'EKG'});
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
        
        % make epoch with all R-peaks
        EEG = pop_epoch( EEG, {  'Rpeak'  }, [-0.15   1.5], 'newname', 'Rpeak epochs all', 'epochinfo', 'yes');
        EEG = pop_rmbase( EEG, [-150 -50] ,[]);
             
%         % in order to ensure that each epoch is associated with one
%         % T-locked HEP and one IBI we need to remove teh second T event
%         % from each R-locked epoch
%         Tevent = []; 
%         for evnt = 1:length(EEG.event)
%             if strcmp(EEG.event(evnt).type,  'Tpeak')
%                 Tevent = [Tevent; evnt EEG.event(evnt).epoch];
%             end
%         end
%         evnt2del =[];
%         for evnt = 2:size(Tevent, 1)
%              if Tevent(evnt, 2) == Tevent(evnt-1, 2)
%                     evnt2del = [evnt2del; Tevent(evnt, 1)];
%              end
%         end
%         EEG.event(evnt2del) = [];
        
        t_peaks = 0;
        for evn = 1:length(EEG.event)
            if strcmp(EEG.event(evn).type, 'Tpeak')
                t_peaks = t_peaks + 1;
            end
        end
        if t_peaks ~= 0
            EEG = pop_epoch( EEG, {  'Tpeak'  }, [-0.3         0.4], 'newname', 'Tpeak epochs all', 'epochinfo', 'yes');
            [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
            if group == 1
                HEP_T_locked_all{group, 2}(yng, t, :, :) = mean(EEG.data, 3);
            else
                HEP_T_locked_all{group, 2}(old, t, :, :) = mean(EEG.data, 3);
            end
        else
            without_t = [without_t; subj_number];
            if group == 1
                HEP_T_locked_all{group, 2}(yng, t, :, :) = NaN;
            else
                HEP_T_locked_all{group, 2}(old, t, :, :) = NaN;
            end
        end
        
%         EEG = pop_delset( EEG, length(ALLEEG) );
        
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
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 3,'gui','off'); 
        EEG = pop_epoch( EEG, {  'Rpeak'  }, [-0.15   .850], 'newname', 'Rpeak epochs b4cue', 'epochinfo', 'yes');
        EEG = pop_rmbase( EEG, [-150 -50] ,[]);
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 3,'gui','off'); 
%         % in order to ensure that each epoch is associated with one
%         % T-locked HEP and one IBI we need to remove teh second T event
%         % from each R-locked epoch
%         Tevent = []; 
%         for evnt = 1:length(EEG.event)
%             if strcmp(EEG.event(evnt).type,  'Tpeak')
%                 Tevent = [Tevent; evnt EEG.event(evnt).epoch];
%             end
%         end
%         evnt2del =[];
%         for evnt = 2:size(Tevent, 1)
%              if Tevent(evnt, 2) == Tevent(evnt-1, 2)
%                     evnt2del = [evnt2del; Tevent(evnt, 1)];
%              end
%         end
%         EEG.event(evnt2del) = [];
        
        t_peaks = 0;
        for evn = 1:length(EEG.event)
            if strcmp(EEG.event(evn).type, 'Tpeak')
                t_peaks = t_peaks + 1;
            end
        end
        if t_peaks ~= 0 % in case no t peaks were detected
            EEG = pop_epoch( EEG, {  'Tpeak'  }, [-0.3         0.4], 'newname', 'Tpeak epochs b4cue', 'epochinfo', 'yes');
%             [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'gui','off'); 
            number_epochs_b4cue(s, t+1) = size(EEG.data, 3);
            if group == 1
                HEP_T_locked_b4cue{group, 2}(yng, t, :, :) = mean(EEG.data, 3);
            else
                HEP_T_locked_b4cue{group, 2}(old, t, :, :) = mean(EEG.data, 3);
            end
        else
            if group == 1
                HEP_T_locked_b4cue{group, 2}(yng, t, :, :) = NaN;
            else
                HEP_T_locked_b4cue{group, 2}(old, t, :, :) = NaN;
            end
        end

%         EEG = pop_delset( EEG, [length(ALLEEG)-1, length(ALLEEG)] );
        
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
        
        EEG = pop_epoch( EEG, {  '1'  }, [.350     2], 'newname', 'epochs aftercue', 'epochinfo', 'yes');
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, length(ALLEEG)+1,'gui','off'); 
        EEG = pop_epoch( EEG, {  'Rpeak'  }, [-0.15     .850], 'newname', 'Rpeak epochs aftercue', 'epochinfo', 'yes');
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, length(ALLEEG)+1,'gui','off'); 
        EEG = pop_rmbase( EEG, [-150 -50] , []);
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 3,'gui','off'); 
        
%         % in order to ensure that each epoch is associated with one
%         % T-locked HEP and one IBI we need to remove teh second T event
%         % from each R-locked epoch
%         Tevent = []; 
%         for evnt = 1:length(EEG.event)
%             if strcmp(EEG.event(evnt).type,  'Tpeak')
%                 Tevent = [Tevent; evnt EEG.event(evnt).epoch];
%             end
%         end
%         evnt2del =[];
%         for evnt = 2:size(Tevent, 1)
%              if Tevent(evnt, 2) == Tevent(evnt-1, 2)
%                     evnt2del = [evnt2del; Tevent(evnt, 1)];
%              end
%         end
%         EEG.event(evnt2del) = [];
        
        t_peaks = 0;
        for evn = 1:length(EEG.event)
            if strcmp(EEG.event(evn).type, 'Tpeak')
                t_peaks = t_peaks + 1;
            end
        end
        if t_peaks ~= 0
            EEG = pop_epoch( EEG, {  'Tpeak'  }, [-0.3  0.4], 'newname', 'Tpeak epochs aftercue', 'epochinfo', 'yes');
            [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
            number_epochs_after_cue(s, t+1) = size(EEG.data, 3);
            if group == 1
                HEP_T_locked_aftercue{group, 2}(yng, t, :, :) = mean(EEG.data, 3);
            else
                HEP_T_locked_aftercue{group, 2}(old, t, :, :) = mean(EEG.data, 3);
            end
        else
            if group == 1
                HEP_T_locked_aftercue{group, 2}(yng, t, :, :) = NaN;
            else
                HEP_T_locked_aftercue{group, 2}(old, t, :, :) = NaN;
            end
        end
    end       
end

save HEP_T_locked_aftercue HEP_T_locked_aftercue
save HEP_T_locked_b4cue HEP_T_locked_b4cue
save HEP_T_locked_all HEP_T_locked_all

save number_epochs_after_cue number_epochs_after_cue
save number_epochs_b4cue number_epochs_b4cue

%% plot figures

channels = load('chanlocs.mat');
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

        x_axis=-0.299:0.002:0.4; % in seconds
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
        if t == 1
            axis([-inf inf -2 2])
        else
            axis([-inf inf -3 1])
        end
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
        if t == 1
            axis([-inf inf -2 2])
        else
            axis([-inf inf -3 1])
        end
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
        if t == 1
            axis([-inf inf -2 2])
        else
            axis([-inf inf -3 1])
        end
        ax.LineWidth = 2.5;
        ax.FontSize = 28;
        ax.FontName = 'Arial';
        xlabel('Time (s)', 'FontSize', 32, 'FontWeight','normal')
        ylabel('Amplitude (\muV)', 'FontSize', 32, 'FontWeight','normal')
        name_channels =  channels.chanlocs(p).labels;
        title(['HEP_T_locked ', tasks{t}, name_channels] , 'FontSize', 32, 'FontWeight','normal')

    end

end

%% mean and sd of number of epochs used
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


%%
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

