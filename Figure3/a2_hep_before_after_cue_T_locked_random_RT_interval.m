% calculate HEP_T_locked 2s interval before cue and 2 s interval after cue
clear; close all;
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
older = [11 12 14 20 21 22 23 32 37 38 41 43 47 48 49 52 55 57 58 63 64 65 67 69 7 70 71 75 8 83 86];
young = [13 15 16 25 26 28 31 33 34 36 4 42 44 45 46 50 51 53 54 56 59 6 62 66 68 72 74 76 78 80 82 84 85 9];

% note: AB21, AB50 and AB55 were excluded from T peak analyses because (s = 8 30 35)
% these were not clearly visible in the ECGs

% load RT intervals for each participant
addpath('E:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\LIMO_stats\HEP_grp_task_effect_Tlocked');
% RT = readtable("E:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\LIMO_stats\HEP_grp_task_effect_Tlocked\RT_interval_2.xls");
% RT = readtable('E:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\LIMO_stats\HEP_grp_task_effect_Tlocked\RT_interval_2.csv');
RT_dir = 'E:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\LIMO_stats\HEP_grp_task_effect_Tlocked';
load([RT_dir filesep 'RT_interval']);
% all files
files_folder = 'E:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\eeg_derivatives\random_R_T_events';
folder_content = dir(files_folder);

% subjects included
subjects = {};
for f = 1:length(folder_content)
    if contains(folder_content(f).name, 'eeg_random_Rpeak_Tpeak')
        subjects = [subjects, folder_content(f).name(28:31)];
    end
end
subjects = unique(subjects);
% subjects excluded from T-locked HEP analyses due to difficulty
% identifying AB21, AB50 and AB55
subjects = setdiff(subjects, {'AB21', 'AB50', 'AB55'});

tasks = {'passive', 'simpleRT', 'gonogo'};

HEP_T_locked_random = cell(2, 2); HEP_T_locked_aftercue_random = cell(2, 2); HEP_T_locked_b4cue_random = cell(2,2);
HEP_R_locked_random = cell(2, 2);

number_trials_random_T = cell(2); % group - number of trials: all, after cue, b4 cue
number_trials_random_T{1} = zeros(length(young), 3);
number_trials_random_T{2} = zeros(length(older), 3);

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
    
    HEP_T_locked_random{group, 1} = [HEP_T_locked_random{group, 1} subj_number];
    HEP_T_locked_aftercue_random{group, 1} = [HEP_T_locked_aftercue_random{group, 1}; subj_number];
    HEP_T_locked_b4cue_random{group, 1} = [HEP_T_locked_b4cue_random{group, 1}; subj_number];
    HEP_R_locked_random{group, 1} = [HEP_R_locked_random{group, 1}; subj_number];
       
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
        EEG = pop_rmbase( EEG, [-150 -50] ,[]);
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
        % save data with HEP locked with random R events
        if group == 1
            HEP_R_locked_random{group, 2}(yng, t, :, :) = mean(EEG.data, 3);
        else
            HEP_R_locked_random{group, 2}(old, t, :, :) = mean(EEG.data, 3);
        end
        
        
        RT_intrvl = [];
        subj_idx = find(RT_interval{group, 1}(:, 1) == subj_number);
        if t == 1
            RT_intrvl = mean([RT_interval{group, 1}(subj_idx, 2), RT_interval{group, 2}(subj_idx, 2)]);
%             RT_intrvl = mean([RT.RT_interv_simpleRT(RT.participant == subj_number), RT.RT_interv_gng(RT.participant == subj_number)]);
        elseif t == 2
            RT_intrvl = RT_interval{group, 1}(subj_idx, 2);
        elseif t == 3
            RT_intrvl = RT_interval{group, 2}(subj_idx, 2);
        end
        % cut epoch locked to R random plus the avg RT interval of that
        % participant
        EEG = pop_select( EEG, 'time',[(RT_intrvl-300)*.001  (RT_intrvl+400)*.001]);
%         EEG = pop_epoch( EEG, {  'Tpeak'  }, [-0.2         0.4], 'newname', 'Tpeak epochs all', 'epochinfo', 'yes');
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
        if group == 1
            number_trials_random_T{1}(yng, 1) = size(EEG.data, 3);
            HEP_T_locked_random{group, 2}(yng, t, :, :) = mean(EEG.data, 3);
        else
            number_trials_random_T{2}(old, 1) = size(EEG.data, 3);
            HEP_T_locked_random{group, 2}(old, t, :, :) = mean(EEG.data, 3);
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
        EEG = pop_epoch( EEG, {  '1'  }, [-1.65        0], 'newname', 'epochs b4cue', 'epochinfo', 'yes');
        EEG = pop_epoch( EEG, {  'Rpeak'  }, [-0.15         .85], 'newname', 'Rpeak epochs b4cue', 'epochinfo', 'yes');
        EEG = pop_rmbase( EEG, [-150 -50] ,[]);
        
        EEG = pop_select( EEG, 'time',[(RT_intrvl-300)*.001  (RT_intrvl+400)*.001]);
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
%         [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'gui','off'); 
        if group == 1
            number_trials_random_T{1}(yng, 2) = size(EEG.data, 3);
            HEP_T_locked_b4cue_random{group, 2}(yng, t, :, :) = mean(EEG.data, 3);
        else
            number_trials_random_T{2}(old, 2) = size(EEG.data, 3);
            HEP_T_locked_b4cue_random{group, 2}(old, t, :, :) = mean(EEG.data, 3);
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
        EEG = pop_epoch( EEG, {  '1'  }, [0     2], 'newname', 'epochs aftercue', 'epochinfo', 'yes');
        EEG = pop_select( EEG, 'time',[.350 2]);
        EEG = pop_epoch( EEG, {  'Rpeak'  }, [-0.15     .85], 'newname', 'Rpeak epochs aftercue', 'epochinfo', 'yes');
        EEG = pop_rmbase( EEG, [-150 -50] ,[]);
        
        EEG = pop_select( EEG, 'time',[(RT_intrvl-300)*.001  (RT_intrvl+400)*.001]);
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
        if group == 1
            number_trials_random_T{1}(yng, 3) = size(EEG.data, 3);
            HEP_T_locked_aftercue_random{group, 2}(yng, t, :, :) = mean(EEG.data, 3);
        else
            number_trials_random_T{2}(old, 3) = size(EEG.data, 3);
            HEP_T_locked_aftercue_random{group, 2}(old, t, :, :) = mean(EEG.data, 3);
        end

    end       
end

save_dir = 'E:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\results_random_events';
cd(save_dir)
save HEP_T_locked_aftercue_random HEP_T_locked_aftercue_random
save HEP_T_locked_b4cue_random HEP_T_locked_b4cue_random
save HEP_T_locked_random HEP_T_locked_random
save number_trials_random_T number_trials_random_T
save HEP_R_locked_random HEP_R_locked_random

%% plot figures
% R-HEP from random R events - all cycles
save_dir = 'E:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\results_random_events';
cd(save_dir)
chan_dir = 'E:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study';
channels = load([chan_dir filesep 'chanlocs.mat']);
tasks = {'passive', 'simpleRT', 'gonogo'};
load HEP_R_locked_random;

x_axis = -0.149:.002:1.5;

for t = 1:length(tasks)
    % HEP_R_locked all
    young_all_mean = squeeze(mean(HEP_R_locked_random{1, 2}(:, t, :, :), 1, 'omitnan'));
    older_all_mean = squeeze(mean(HEP_R_locked_random{2, 2}(:, t, :, :), 1, 'omitnan'));

    young_all_se = squeeze(std(HEP_R_locked_random{1, 2}(:, t, :, :), [], 1, 'omitnan'))/sqrt(size(HEP_R_locked_random{1, 2}, 1));
    older_all_se = squeeze(std(HEP_R_locked_random{2, 2}(:, t, :, :), [], 1, 'omitnan'))/sqrt(size(HEP_R_locked_random{1, 2}, 1));
    channels_num = 51%34;%[7, 16, 25, 34];
    % close all
    for p = channels_num
        figure;
        plot(x_axis, young_all_mean(p, :), 'color', 'k',  'LineWidth',1.5)
        hold on
        jbfill(x_axis, [young_all_mean(p, :) + young_all_se(p, :)], [young_all_mean(p, :) - young_all_se(p, :)], 'k', 'k', 1, .2)
        hold on
        plot(x_axis, older_all_mean(p, :), '--r',  'LineWidth',1.5) 
        hold on
        jbfill(x_axis, [older_all_mean(p, :) + older_all_se(p, :)], [older_all_mean(p, :) - older_all_se(p, :)], 'r', 'r', 1, .2)
        hold on
        plot(x_axis, zeros(length(x_axis), 1), ':k',  'LineWidth',1) 
        hold off
        box off
        ax = gca;
        if t == 1
            axis([-inf inf -2 2])
        else
            axis([.2 .6 -1 1])
        end
        ax.LineWidth = 2.5;
        ax.FontSize = 28;
        ax.FontName = 'Arial';
        xlabel('Time (s)', 'FontSize', 32, 'FontWeight','normal')
        ylabel('Amplitude (\muV)', 'FontSize', 32, 'FontWeight','normal')
        name_channels =  channels.chanlocs(p).labels;
        title(['R-HEP random ', tasks{t}, name_channels] , 'FontSize', 32, 'FontWeight','normal')

    end
end

%% T-HEP from random R events - all cycles
save_dir = 'E:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\results_random_events';
cd(save_dir)
chan_dir = 'E:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study';
channels = load([chan_dir filesep 'chanlocs.mat']);
tasks = {'passive', 'simpleRT', 'gonogo'};
load HEP_T_locked_random;

x_axis = -.299:.002:.4;

for t = 1:length(tasks)
    % HEP_R_locked all
    young_all_mean = squeeze(mean(HEP_T_locked_random{1, 2}(:, t, :, :), 1, 'omitnan'));
    older_all_mean = squeeze(mean(HEP_T_locked_random{2, 2}(:, t, :, :), 1, 'omitnan'));

    young_all_se = squeeze(std(HEP_T_locked_random{1, 2}(:, t, :, :), [], 1, 'omitnan'))/sqrt(size(HEP_T_locked_random{1, 2}, 1));
    older_all_se = squeeze(std(HEP_T_locked_random{2, 2}(:, t, :, :), [], 1, 'omitnan'))/sqrt(size(HEP_T_locked_random{1, 2}, 1));
    channels_num = 34;%[7, 16, 25, 34];
    % close all
    for p = channels_num
        figure;
        plot(x_axis, young_all_mean(p, :), 'color', 'k',  'LineWidth',1.5)
        hold on
        jbfill(x_axis, [young_all_mean(p, :) + young_all_se(p, :)], [young_all_mean(p, :) - young_all_se(p, :)], 'k', 'k', 1, .2)
        hold on
        plot(x_axis, older_all_mean(p, :), '--r',  'LineWidth',1.5) 
        hold on
        jbfill(x_axis, [older_all_mean(p, :) + older_all_se(p, :)], [older_all_mean(p, :) - older_all_se(p, :)], 'r', 'r', 1, .2)
        hold on
        plot(x_axis, zeros(length(x_axis), 1), ':k',  'LineWidth',1) 
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
        title(['T-HEP random ', tasks{t}, name_channels] , 'FontSize', 32, 'FontWeight','normal')

    end
end

%% T-HEP from random R events - all cycles - both groups and tasks together
save_dir = 'E:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\results_random_events';
cd(save_dir)
chan_dir = 'E:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study';
channels = load([chan_dir filesep 'chanlocs.mat']);
tasks = {'passive', 'simpleRT', 'gonogo'};
load HEP_T_locked_random;
data = cat(1, HEP_T_locked_random{1, 2}, HEP_T_locked_random{2, 2}); % participants, task, channels, time
data = mean(data(:, 2:3, :, :), 2);

data_all_mean = squeeze(mean(data, 1, 'omitnan'));

data_all_se = squeeze(std(data, [], 1, 'omitnan'))/sqrt(size(data, 1));
x_axis = -299:2:400;
channels_num = [7, 16, 25, 34, 51];
% close all
for p = channels_num
    figure;
    plot(x_axis, data_all_mean(p, :), 'color', 'k',  'LineWidth',1.5)
    hold on
    jbfill(x_axis, [data_all_mean(p, :) + data_all_se(p, :)], [data_all_mean(p, :) - data_all_se(p, :)], 'k', 'k', 1, .2)
    hold on
    plot(x_axis, zeros(length(x_axis), 1), ':k',  'LineWidth',1) 
    hold off
    box off
    ax = gca;
    axis([-inf inf -1.5 2.25])
    ax.LineWidth = 2.5;
    ax.FontSize = 28;
    ax.FontName = 'Arial';
    xlabel('Time (ms)', 'FontSize', 32, 'FontWeight','normal')
    ylabel('Amplitude (\muV)', 'FontSize', 32, 'FontWeight','normal')
    name_channels =  channels.chanlocs(p).labels;
    title(name_channels , 'FontSize', 32, 'FontWeight','normal')

end


%%
save_dir = 'E:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\results_random_events';
cd(save_dir)
chan_dir = 'E:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study';
channels = load([chan_dir filesep 'chanlocs.mat']);
tasks = {'passive', 'simple RT', 'go/no-go'};
load HEP_T_locked_aftercue_random; load HEP_T_locked_b4cue_random; load HEP_T_locked_random;

for t = 1:length(tasks)
% 
% figure;
% for p = 1:size(HEP_T_locked_b4cue{2, 2}, 1)
%     plot(squeeze(HEP_T_locked_b4cue{2, 2}(p, t, 7, :)))
%     waitforbuttonpress
% end
% close all

% HEP_T_locked before cue
young_b4cue_mean = squeeze(mean(HEP_T_locked_b4cue_random{1, 2}(:, t, :, :), 1, 'omitnan'));
older_b4cue_mean = squeeze(mean(HEP_T_locked_b4cue_random{2, 2}(:, t, :, :), 1, 'omitnan'));

young_b4cue_se = squeeze(std(HEP_T_locked_b4cue_random{1, 2}(:, t, :, :), [], 1, 'omitnan'))/sqrt(size(HEP_T_locked_b4cue_random{1, 2}, 1));
older_b4cue_se = squeeze(std(HEP_T_locked_b4cue_random{2, 2}(:, t, :, :), [], 1, 'omitnan'))/sqrt(size(HEP_T_locked_b4cue_random{1, 2}, 1));

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
young_aftercue_mean = squeeze(mean(HEP_T_locked_aftercue_random{1, 2}(:, t, :, :), 1, 'omitnan'));
older_aftercue_mean = squeeze(mean(HEP_T_locked_aftercue_random{2, 2}(:, t, :, :), 1, 'omitnan'));

young_aftercue_se = squeeze(std(HEP_T_locked_aftercue_random{1, 2}(:, t, :, :), [], 1, 'omitnan'))/sqrt(size(HEP_T_locked_aftercue_random{1, 2}, 1));
older_aftercue_se = squeeze(std(HEP_T_locked_aftercue_random{2, 2}(:, t, :, :), [], 1, 'omitnan'))/sqrt(size(HEP_T_locked_aftercue_random{1, 2}, 1));

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
young_all_mean = squeeze(mean(HEP_T_locked_random{1, 2}(:, t, :, :), 1, 'omitnan'));
older_all_mean = squeeze(mean(HEP_T_locked_random{2, 2}(:, t, :, :), 1, 'omitnan'));

young_all_se = squeeze(std(HEP_T_locked_random{1, 2}(:, t, :, :), [], 1, 'omitnan'))/sqrt(size(HEP_T_locked_random{1, 2}, 1));
older_all_se = squeeze(std(HEP_T_locked_random{2, 2}(:, t, :, :), [], 1, 'omitnan'))/sqrt(size(HEP_T_locked_random{1, 2}, 1));


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


channels_num = 34;%[7, 16, 25, 34];
% close all
    for p = channels_num

        x_axis=-0.299:0.002:0.4; % in seconds
        figure;
        plot(x_axis, young_b4cue_mean(p, :), '--', 'color', [0 0.4470 0.7410],  'LineWidth',1.5) %  (light blue)
        hold on
        jbfill(x_axis, [young_b4cue_mean(p, :) + young_b4cue_se(p, :)], [young_b4cue_mean(p, :) - young_b4cue_se(p, :)], [0 0.4470 0.7410],[0 0.4470 0.7410], 1, 0.2)
        hold on
        plot(x_axis, young_aftercue_mean(p, :), 'color', 'k',  'LineWidth',1.5)
        hold on
        jbfill(x_axis, [young_aftercue_mean(p, :) + young_aftercue_se(p, :)], [young_aftercue_mean(p, :) - young_aftercue_se(p, :)], 'k', 'k', 1, 0.2)
        hold on
        plot(x_axis, zeros(length(x_axis), 1), ':k')
        hold off
        box off
        ax = gca;
        if t == 1
            axis([-inf inf -2 2])
        else
            axis([-inf inf -5 1])
        end
        ax.LineWidth = 2.5;
        ax.FontSize = 28;
        ax.FontName = 'Arial';
        xlabel('Time (s)', 'FontSize', 32, 'FontWeight','normal')
        ylabel('Amplitude (\muV)', 'FontSize', 32, 'FontWeight','normal')
        name_channels =  channels.chanlocs(p).labels;
        title(['Young - ', tasks{t}, ' - ', name_channels] , 'FontSize', 32, 'FontWeight','normal')
%         basefilename2 = sprintf('Tepochplot_gng-sRT_young_cut_%d.fig', p);
%         savefig(basefilename2);


        figure;
        plot(x_axis, older_b4cue_mean(p, :), '--', 'color', [0 0.4470 0.7410],  'LineWidth',1.5) %plot mean gng (light blue)
        hold on
        jbfill(x_axis, [older_b4cue_mean(p, :) + older_b4cue_se(p, :)], [older_b4cue_mean(p, :) - older_b4cue_se(p, :)], [0 0.4470 0.7410],[0 0.4470 0.7410], 1, 0.2)
        hold on
        plot(x_axis, older_aftercue_mean(p, :), 'color', 'k',  'LineWidth',1.5) %plot mean gng (light blue)
        hold on
        jbfill(x_axis, [older_aftercue_mean(p, :) + older_aftercue_se(p, :)], [older_aftercue_mean(p, :) - older_aftercue_se(p, :)], 'k', 'k', 1, 0.2)
        hold on
        plot(x_axis, zeros(length(x_axis), 1), ':k')
        hold off
        box off
        ax = gca;
        if t == 1
            axis([-inf inf -2 2])
        else
            axis([-inf inf -5 1])
        end
        ax.LineWidth = 2.5;
        ax.FontSize = 28;
        ax.FontName = 'Arial';
        xlabel('Time (s)', 'FontSize', 32, 'FontWeight','normal')
        ylabel('Amplitude (\muV)', 'FontSize', 32, 'FontWeight','normal')
        name_channels =  channels.chanlocs(p).labels;
        title(['Older - ', tasks{t}, ' - ', name_channels] , 'FontSize', 32, 'FontWeight','normal')
        
        
        figure;
        plot(x_axis, young_all_mean(p, :), 'color', 'k',  'LineWidth',1.5)
        hold on
        jbfill(x_axis, [young_all_mean(p, :) + young_all_se(p, :)], [young_all_mean(p, :) - young_all_se(p, :)], 'k', 'k')
        hold on
        plot(x_axis, older_all_mean(p, :), '--r',  'LineWidth',1.5) 
        hold on
        jbfill(x_axis, [older_all_mean(p, :) + older_all_se(p, :)], [older_all_mean(p, :) - older_all_se(p, :)], 'r', 'r')
        hold on
        plot(x_axis, zeros(length(x_axis), 1), ':k')
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
        title(['T-HEP random ', tasks{t}, name_channels] , 'FontSize', 32, 'FontWeight','normal')

    end
end

%% plot T-HEP random with average on both task conditions 
save_dir = 'E:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\results_random_events';
cd(save_dir)
chan_dir = 'E:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study';
channels = load([chan_dir filesep 'chanlocs.mat']);
tasks = {'passive', 'simple RT', 'go/no-go'};
load HEP_T_locked_aftercue_random; load HEP_T_locked_b4cue_random; load HEP_T_locked_random;


% HEP_T_locked before cue
young_b4cue_mean = squeeze(mean(mean(HEP_T_locked_b4cue_random{1, 2}(:, 2:3, :, :), 2, 'omitnan'), 1, 'omitnan'));
older_b4cue_mean = squeeze(mean(mean(HEP_T_locked_b4cue_random{2, 2}(:, 2:3, :, :), 2, 'omitnan'), 1, 'omitnan'));

young_b4cue_se = squeeze(std(mean(HEP_T_locked_b4cue_random{1, 2}(:, 2:3, :, :), 2, 'omitnan'), [], 1, 'omitnan'))/sqrt(size(HEP_T_locked_b4cue_random{1, 2}, 1));
older_b4cue_se = squeeze(std(mean(HEP_T_locked_b4cue_random{2, 2}(:, 2:3, :, :), 2, 'omitnan'), [], 1, 'omitnan'))/sqrt(size(HEP_T_locked_b4cue_random{1, 2}, 1));


% HEP_T_locked after cue
young_aftercue_mean = squeeze(mean(mean(HEP_T_locked_aftercue_random{1, 2}(:, 2:3, :, :), 2, 'omitnan'), 1, 'omitnan'));
older_aftercue_mean = squeeze(mean(mean(HEP_T_locked_aftercue_random{2, 2}(:, 2:3, :, :), 2, 'omitnan'), 1, 'omitnan'));

young_aftercue_se = squeeze(std(mean(HEP_T_locked_aftercue_random{1, 2}(:, 2:3, :, :), 2, 'omitnan'), [], 1, 'omitnan'))/sqrt(size(HEP_T_locked_aftercue_random{1, 2}, 1));
older_aftercue_se = squeeze(std(mean(HEP_T_locked_aftercue_random{2, 2}(:, 2:3, :, :), 2, 'omitnan'), [], 1, 'omitnan'))/sqrt(size(HEP_T_locked_aftercue_random{1, 2}, 1));


channels_num = 34;%[7, 16, 25, 34];
% close all
    for p = channels_num

        x_axis=-0.299:0.002:0.4; % in seconds
        figure;
        plot(x_axis, young_b4cue_mean(p, :), '--', 'color', [0 0.4470 0.7410],  'LineWidth',1.5) %  (light blue)
        hold on
        jbfill(x_axis, [young_b4cue_mean(p, :) + young_b4cue_se(p, :)], [young_b4cue_mean(p, :) - young_b4cue_se(p, :)], [0 0.4470 0.7410],[0 0.4470 0.7410], 1, 0.2)
        hold on
        plot(x_axis, young_aftercue_mean(p, :), 'color', 'k',  'LineWidth',1.5)
        hold on
        jbfill(x_axis, [young_aftercue_mean(p, :) + young_aftercue_se(p, :)], [young_aftercue_mean(p, :) - young_aftercue_se(p, :)], 'k', 'k', 1, 0.2)
        hold on
        plot(x_axis, zeros(length(x_axis), 1), ':k')
        hold off
        box off
        ax = gca;
        if t == 1
            axis([-inf inf -2 2])
        else
            axis([-inf inf -5 1])
        end
        ax.LineWidth = 2.5;
        ax.FontSize = 28;
        ax.FontName = 'Arial';
        xlabel('Time (s)', 'FontSize', 32, 'FontWeight','normal')
        ylabel('Amplitude (\muV)', 'FontSize', 32, 'FontWeight','normal')
        name_channels =  channels.chanlocs(p).labels;
        title(['Young - ', name_channels] , 'FontSize', 32, 'FontWeight','normal')
%         basefilename2 = sprintf('Tepochplot_gng-sRT_young_cut_%d.fig', p);
%         savefig(basefilename2);


        figure;
        plot(x_axis, older_b4cue_mean(p, :), '--', 'color', [0 0.4470 0.7410],  'LineWidth',1.5) %plot mean gng (light blue)
        hold on
        jbfill(x_axis, [older_b4cue_mean(p, :) + older_b4cue_se(p, :)], [older_b4cue_mean(p, :) - older_b4cue_se(p, :)], [0 0.4470 0.7410],[0 0.4470 0.7410], 1, 0.2)
        hold on
        plot(x_axis, older_aftercue_mean(p, :), 'color', 'k',  'LineWidth',1.5) %plot mean gng (light blue)
        hold on
        jbfill(x_axis, [older_aftercue_mean(p, :) + older_aftercue_se(p, :)], [older_aftercue_mean(p, :) - older_aftercue_se(p, :)], 'k', 'k', 1, 0.2)
        hold on
        plot(x_axis, zeros(length(x_axis), 1), ':k')
        hold off
        box off
        ax = gca;
        if t == 1
            axis([-inf inf -2 2])
        else
            axis([-inf inf -5 1])
        end
        ax.LineWidth = 2.5;
        ax.FontSize = 28;
        ax.FontName = 'Arial';
        xlabel('Time (s)', 'FontSize', 32, 'FontWeight','normal')
        ylabel('Amplitude (\muV)', 'FontSize', 32, 'FontWeight','normal')
        name_channels =  channels.chanlocs(p).labels;
        title(['Older - ', name_channels] , 'FontSize', 32, 'FontWeight','normal')

    end


    
 %% plot T-HEP b4 and after cue with average on both task conditions 

chan_dir = 'E:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study';
cd(chan_dir)
channels = load([chan_dir filesep 'chanlocs.mat']);
tasks = {'passive', 'simple RT', 'go/no-go'};
load HEP_T_locked_aftercue; load HEP_T_locked_b4cue;


% HEP_T_locked before cue
young_b4cue_mean = squeeze(mean(mean(HEP_T_locked_b4cue{1, 2}(:, 2:3, :, :), 2, 'omitnan'), 1, 'omitnan'));
older_b4cue_mean = squeeze(mean(mean(HEP_T_locked_b4cue{2, 2}(:, 2:3, :, :), 2, 'omitnan'), 1, 'omitnan'));

young_b4cue_se = squeeze(std(mean(HEP_T_locked_b4cue{1, 2}(:, 2:3, :, :), 2, 'omitnan'), [], 1, 'omitnan'))/sqrt(size(HEP_T_locked_b4cue{1, 2}, 1));
older_b4cue_se = squeeze(std(mean(HEP_T_locked_b4cue{2, 2}(:, 2:3, :, :), 2, 'omitnan'), [], 1, 'omitnan'))/sqrt(size(HEP_T_locked_b4cue{1, 2}, 1));


% HEP_T_locked after cue
young_aftercue_mean = squeeze(mean(mean(HEP_T_locked_aftercue{1, 2}(:, 2:3, :, :), 2, 'omitnan'), 1, 'omitnan'));
older_aftercue_mean = squeeze(mean(mean(HEP_T_locked_aftercue{2, 2}(:, 2:3, :, :), 2, 'omitnan'), 1, 'omitnan'));

young_aftercue_se = squeeze(std(mean(HEP_T_locked_aftercue{1, 2}(:, 2:3, :, :), 2, 'omitnan'), [], 1, 'omitnan'))/sqrt(size(HEP_T_locked_aftercue{1, 2}, 1));
older_aftercue_se = squeeze(std(mean(HEP_T_locked_aftercue{2, 2}(:, 2:3, :, :), 2, 'omitnan'), [], 1, 'omitnan'))/sqrt(size(HEP_T_locked_aftercue{1, 2}, 1));


channels_num = 34;%[7, 16, 25, 34];
% close all
    for p = channels_num

        x_axis=-0.299:0.002:0.4; % in seconds
        figure;
        plot(x_axis, young_b4cue_mean(p, :), '--', 'color', [0 0.4470 0.7410],  'LineWidth',1.5) %  (light blue)
        hold on
        jbfill(x_axis, [young_b4cue_mean(p, :) + young_b4cue_se(p, :)], [young_b4cue_mean(p, :) - young_b4cue_se(p, :)], [0 0.4470 0.7410],[0 0.4470 0.7410], 1, 0.2)
        hold on
        plot(x_axis, young_aftercue_mean(p, :), 'color', 'k',  'LineWidth',1.5)
        hold on
        jbfill(x_axis, [young_aftercue_mean(p, :) + young_aftercue_se(p, :)], [young_aftercue_mean(p, :) - young_aftercue_se(p, :)], 'k', 'k', 1, 0.2)
        hold on
        plot(x_axis, zeros(length(x_axis), 1), ':k')
        hold off
        box off
        ax = gca;
        if t == 1
            axis([-inf inf -2 2])
        else
            axis([-inf inf -5 1])
        end
        ax.LineWidth = 2.5;
        ax.FontSize = 28;
        ax.FontName = 'Arial';
        xlabel('Time (s)', 'FontSize', 32, 'FontWeight','normal')
        ylabel('Amplitude (\muV)', 'FontSize', 32, 'FontWeight','normal')
        name_channels =  channels.chanlocs(p).labels;
        title(['Young - ', name_channels] , 'FontSize', 32, 'FontWeight','normal')
%         basefilename2 = sprintf('Tepochplot_gng-sRT_young_cut_%d.fig', p);
%         savefig(basefilename2);


        figure;
        plot(x_axis, older_b4cue_mean(p, :), '--', 'color', [0 0.4470 0.7410],  'LineWidth',1.5) %plot mean gng (light blue)
        hold on
        jbfill(x_axis, [older_b4cue_mean(p, :) + older_b4cue_se(p, :)], [older_b4cue_mean(p, :) - older_b4cue_se(p, :)], [0 0.4470 0.7410],[0 0.4470 0.7410], 1, 0.2)
        hold on
        plot(x_axis, older_aftercue_mean(p, :), 'color', 'k',  'LineWidth',1.5) %plot mean gng (light blue)
        hold on
        jbfill(x_axis, [older_aftercue_mean(p, :) + older_aftercue_se(p, :)], [older_aftercue_mean(p, :) - older_aftercue_se(p, :)], 'k', 'k', 1, 0.2)
        hold on
        plot(x_axis, zeros(length(x_axis), 1), ':k')
        hold off
        box off
        ax = gca;
        if t == 1
            axis([-inf inf -2 2])
        else
            axis([-inf inf -5 1])
        end
        ax.LineWidth = 2.5;
        ax.FontSize = 28;
        ax.FontName = 'Arial';
        xlabel('Time (s)', 'FontSize', 32, 'FontWeight','normal')
        ylabel('Amplitude (\muV)', 'FontSize', 32, 'FontWeight','normal')
        name_channels =  channels.chanlocs(p).labels;
        title(['Older - ', name_channels] , 'FontSize', 32, 'FontWeight','normal')

    end   
    
    
%% plot T-HEP b4 and after cue with average on both task conditions and average across groups

chan_dir = 'E:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study';
cd(chan_dir)
channels = load([chan_dir filesep 'chanlocs.mat']);
tasks = {'passive', 'simple RT', 'go/no-go'};
load HEP_T_locked_aftercue; load HEP_T_locked_b4cue;


% HEP_T_locked before cue
b4cue_mean = squeeze(mean(mean(cat(1, HEP_T_locked_b4cue{1, 2}(:, 2:3, :, :),...
   HEP_T_locked_b4cue{2, 2}(:, 2:3, :, :)), 2, 'omitnan'), 1, 'omitnan'));

b4cue_se = squeeze(std(mean(cat(1, HEP_T_locked_b4cue{1, 2}(:, 2:3, :, :),...
   HEP_T_locked_b4cue{2, 2}(:, 2:3, :, :)), 2, 'omitnan'), [], 1, 'omitnan'))/sqrt(size(cat(1, HEP_T_locked_b4cue{1, 2}, HEP_T_locked_b4cue{2, 2}), 1));


% HEP_T_locked after cue
aftercue_mean = squeeze(mean(mean(cat(1, HEP_T_locked_aftercue{1, 2}(:, 2:3, :, :), ...
    HEP_T_locked_aftercue{2, 2}(:, 2:3, :, :)), 2, 'omitnan'), 1, 'omitnan'));

aftercue_se = squeeze(std(mean(cat(1, HEP_T_locked_aftercue{1, 2}(:, 2:3, :, :), ...
    HEP_T_locked_aftercue{2, 2}(:, 2:3, :, :)), 2, 'omitnan'), [], 1, 'omitnan'))/sqrt(size(cat(1, HEP_T_locked_aftercue{1, 2}, ...
    HEP_T_locked_aftercue{2, 2}), 1));

channels_num = 34;%[7, 16, 25, 34];
% close all
    for p = channels_num

        x_axis=-0.299:0.002:0.4; % in seconds
        figure;
        plot(x_axis, b4cue_mean(p, :), '--', 'color', [0 0.4470 0.7410],  'LineWidth',1.5) %  (light blue)
        hold on
        jbfill(x_axis, b4cue_mean(p, :) + b4cue_se(p, :), b4cue_mean(p, :) - b4cue_se(p, :), [0 0.4470 0.7410],[0 0.4470 0.7410], 1, 0.2)
        hold on
        plot(x_axis, aftercue_mean(p, :), 'color', 'k',  'LineWidth',1.5)
        hold on
        jbfill(x_axis, aftercue_mean(p, :) + aftercue_se(p, :), aftercue_mean(p, :) - aftercue_se(p, :), 'k', 'k', 1, 0.2)
        hold on
        plot(x_axis, zeros(length(x_axis), 1), ':k')
        hold off
        box off
        ax = gca;
        if t == 1
            axis([-inf inf -2 2])
        else
            axis([-inf inf -5 1])
        end
        ax.LineWidth = 2.5;
        ax.FontSize = 28;
        ax.FontName = 'Arial';
        xlabel('Time (s)', 'FontSize', 32, 'FontWeight','normal')
        ylabel('Amplitude (\muV)', 'FontSize', 32, 'FontWeight','normal')
        name_channels =  channels.chanlocs(p).labels;
        title(name_channels , 'FontSize', 32, 'FontWeight','normal')
%         basefilename2 = sprintf('Tepochplot_gng-sRT_young_cut_%d.fig', p);
%         savefig(basefilename2);

    end

%% smooth data to subtract from HEP

save_dir = 'E:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\results_random_events';
cd(save_dir)
chan_dir = 'E:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study';
channels = load([chan_dir filesep 'chanlocs.mat']);
tasks = {'passive', 'simple RT', 'go/no-go'};
load HEP_T_locked_aftercue_random; load HEP_T_locked_b4cue_random; load HEP_T_locked_random;


% load HEPs
heps_dir ='E:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study';
cd(heps_dir)
load HEP_T_locked_aftercue; load HEP_T_locked_b4cue; load HEP_T_locked_all;

HEP_T_aftercue_detrended = HEP_T_locked_aftercue;
HEP_T_b4cue_detrended = HEP_T_locked_b4cue;

smth_HEP_aftercue_random = HEP_T_locked_aftercue;
smth_HEP_b4cue_random = HEP_T_locked_b4cue;
smth_HEP_random = HEP_T_locked_all;

smth_window_length = 150;

for grp =1:2
    for p = 1:size(HEP_T_locked_b4cue{grp, 2}, 1)
        for t = 1:3
            subject = HEP_T_locked_aftercue_random{grp, 1}(p); % should b the same for HEP_T_locked_b4cue_random
            % smooth ERP obtained with random events after cue
            data_tmp = squeeze(HEP_T_locked_aftercue_random{grp, 2}(p, t, 1:59, :)); % channels x time
            data_tmp_smooth = smoothdata(data_tmp, 2, 'gaussian', smth_window_length);
            smth_HEP_aftercue_random{grp, 2}(p, t, :, :) = data_tmp_smooth;
            
            % subtract smoothed ERP from after cue HEP to correct for
            % baseline shift
            for sbj = 1:size(HEP_T_locked_aftercue{grp, 2}, 1)
                if HEP_T_locked_aftercue{grp, 1}(sbj) == subject
                    data = squeeze(HEP_T_locked_aftercue{grp, 2}(sbj, t, :, :)); % channels x time
                end
            end
            HEP_T_aftercue_detrended{grp, 2}(p, t, :, :) = data - data_tmp_smooth;
            
            % smooth HEP-random before cue
            data_tmp = squeeze(HEP_T_locked_b4cue_random{grp, 2}(p, t, 1:59, :)); % channels x time
            data_tmp_smooth = smoothdata(data_tmp, 2, 'gaussian', smth_window_length);
            smth_HEP_b4cue_random{grp, 2}(p, t, :, :) = data_tmp_smooth;
            
            % correct for baseline shift HEP before cue
            for sbj = 1:size(HEP_T_locked_b4cue{grp, 2}, 1)
                if HEP_T_locked_b4cue{grp, 1}(sbj) == subject
                    data = squeeze(HEP_T_locked_b4cue{grp, 2}(sbj, t, :, :)); % channels x time
                end
            end
            HEP_T_b4cue_detrended{grp, 2}(p, t, :, :) = data - data_tmp_smooth;
            
            data_tmp = squeeze(HEP_T_locked_random{grp, 2}(p, t, 1:59, :)); % channels x time
            data_tmp_smooth = smoothdata(data_tmp, 2, 'gaussian', smth_window_length);
            smth_HEP_random{grp, 2}(p, t, :, :) = data_tmp_smooth;
        end
    end
end

 %% plot smoothed random T-HEP b4 and after cue with average on both task conditions 

chan_dir = 'E:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study';
cd(chan_dir)
channels = load([chan_dir filesep 'chanlocs.mat']);
tasks = {'passive', 'simple RT', 'go/no-go'};



% HEP_T_locked before cue
young_b4cue_mean = squeeze(mean(mean(smth_HEP_b4cue_random{1, 2}(:, 2:3, :, :), 2, 'omitnan'), 1, 'omitnan'));
older_b4cue_mean = squeeze(mean(mean(smth_HEP_b4cue_random{2, 2}(:, 2:3, :, :), 2, 'omitnan'), 1, 'omitnan'));

young_b4cue_se = squeeze(std(mean(smth_HEP_b4cue_random{1, 2}(:, 2:3, :, :), 2, 'omitnan'), [], 1, 'omitnan'))/sqrt(size(smth_HEP_b4cue_random{1, 2}, 1));
older_b4cue_se = squeeze(std(mean(smth_HEP_b4cue_random{2, 2}(:, 2:3, :, :), 2, 'omitnan'), [], 1, 'omitnan'))/sqrt(size(smth_HEP_b4cue_random{1, 2}, 1));


% HEP_T_locked after cue
young_aftercue_mean = squeeze(mean(mean(smth_HEP_aftercue_random{1, 2}(:, 2:3, :, :), 2, 'omitnan'), 1, 'omitnan'));
older_aftercue_mean = squeeze(mean(mean(smth_HEP_aftercue_random{2, 2}(:, 2:3, :, :), 2, 'omitnan'), 1, 'omitnan'));

young_aftercue_se = squeeze(std(mean(smth_HEP_aftercue_random{1, 2}(:, 2:3, :, :), 2, 'omitnan'), [], 1, 'omitnan'))/sqrt(size(smth_HEP_aftercue_random{1, 2}, 1));
older_aftercue_se = squeeze(std(mean(smth_HEP_aftercue_random{2, 2}(:, 2:3, :, :), 2, 'omitnan'), [], 1, 'omitnan'))/sqrt(size(smth_HEP_aftercue_random{1, 2}, 1));


channels_num = 34;%[7, 16, 25, 34];
% close all
for p = channels_num

    x_axis=-0.299:0.002:0.4; % in seconds
    figure;
    plot(x_axis, young_b4cue_mean(p, :), '--', 'color', [0 0.4470 0.7410],  'LineWidth',1.5) %  (light blue)
    hold on
    jbfill(x_axis, [young_b4cue_mean(p, :) + young_b4cue_se(p, :)], [young_b4cue_mean(p, :) - young_b4cue_se(p, :)], [0 0.4470 0.7410],[0 0.4470 0.7410], 1, 0.2)
    hold on
    plot(x_axis, young_aftercue_mean(p, :), 'color', 'k',  'LineWidth',1.5)
    hold on
    jbfill(x_axis, [young_aftercue_mean(p, :) + young_aftercue_se(p, :)], [young_aftercue_mean(p, :) - young_aftercue_se(p, :)], 'k', 'k', 1, 0.2)
    hold on
    plot(x_axis, zeros(length(x_axis), 1), ':k')
    hold off
    box off
    ax = gca;
    if t == 1
        axis([-inf inf -2 2])
    else
        axis([-inf inf -5 1])
    end
    ax.LineWidth = 2.5;
    ax.FontSize = 28;
    ax.FontName = 'Arial';
    xlabel('Time (s)', 'FontSize', 32, 'FontWeight','normal')
    ylabel('Amplitude (\muV)', 'FontSize', 32, 'FontWeight','normal')
    name_channels =  channels.chanlocs(p).labels;
    title(['Young - ', name_channels] , 'FontSize', 32, 'FontWeight','normal')
%         basefilename2 = sprintf('Tepochplot_gng-sRT_young_cut_%d.fig', p);
%         savefig(basefilename2);


    figure;
    plot(x_axis, older_b4cue_mean(p, :), '--', 'color', [0 0.4470 0.7410],  'LineWidth',1.5) %plot mean gng (light blue)
    hold on
    jbfill(x_axis, [older_b4cue_mean(p, :) + older_b4cue_se(p, :)], [older_b4cue_mean(p, :) - older_b4cue_se(p, :)], [0 0.4470 0.7410],[0 0.4470 0.7410], 1, 0.2)
    hold on
    plot(x_axis, older_aftercue_mean(p, :), 'color', 'k',  'LineWidth',1.5) %plot mean gng (light blue)
    hold on
    jbfill(x_axis, [older_aftercue_mean(p, :) + older_aftercue_se(p, :)], [older_aftercue_mean(p, :) - older_aftercue_se(p, :)], 'k', 'k', 1, 0.2)
    hold on
    plot(x_axis, zeros(length(x_axis), 1), ':k')
    hold off
    box off
    ax = gca;
    if t == 1
        axis([-inf inf -2 2])
    else
        axis([-inf inf -5 1])
    end
    ax.LineWidth = 2.5;
    ax.FontSize = 28;
    ax.FontName = 'Arial';
    xlabel('Time (s)', 'FontSize', 32, 'FontWeight','normal')
    ylabel('Amplitude (\muV)', 'FontSize', 32, 'FontWeight','normal')
    name_channels =  channels.chanlocs(p).labels;
    title(['Older - ', name_channels] , 'FontSize', 32, 'FontWeight','normal')

end  


 %% plot detrended T-HEP b4 and after cue with average on both task conditions 

chan_dir = 'E:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study';
cd(chan_dir)
channels = load([chan_dir filesep 'chanlocs.mat']);
tasks = {'passive', 'simple RT', 'go/no-go'};



% HEP_T_locked before cue
young_b4cue_mean = squeeze(mean(mean(HEP_T_b4cue_detrended{1, 2}(:, 2:3, :, :), 2, 'omitnan'), 1, 'omitnan'));
older_b4cue_mean = squeeze(mean(mean(HEP_T_b4cue_detrended{2, 2}(:, 2:3, :, :), 2, 'omitnan'), 1, 'omitnan'));

young_b4cue_se = squeeze(std(mean(HEP_T_b4cue_detrended{1, 2}(:, 2:3, :, :), 2, 'omitnan'), [], 1, 'omitnan'))/sqrt(size(HEP_T_b4cue_detrended{1, 2}, 1));
older_b4cue_se = squeeze(std(mean(HEP_T_b4cue_detrended{2, 2}(:, 2:3, :, :), 2, 'omitnan'), [], 1, 'omitnan'))/sqrt(size(HEP_T_b4cue_detrended{1, 2}, 1));


% HEP_T_locked after cue
young_aftercue_mean = squeeze(mean(mean(HEP_T_aftercue_detrended{1, 2}(:, 2:3, :, :), 2, 'omitnan'), 1, 'omitnan'));
older_aftercue_mean = squeeze(mean(mean(HEP_T_aftercue_detrended{2, 2}(:, 2:3, :, :), 2, 'omitnan'), 1, 'omitnan'));

young_aftercue_se = squeeze(std(mean(HEP_T_aftercue_detrended{1, 2}(:, 2:3, :, :), 2, 'omitnan'), [], 1, 'omitnan'))/sqrt(size(HEP_T_aftercue_detrended{1, 2}, 1));
older_aftercue_se = squeeze(std(mean(HEP_T_aftercue_detrended{2, 2}(:, 2:3, :, :), 2, 'omitnan'), [], 1, 'omitnan'))/sqrt(size(HEP_T_aftercue_detrended{1, 2}, 1));


channels_num = 34;%[7, 16, 25, 34];
% close all
for p = channels_num

    x_axis=-0.299:0.002:0.4; % in seconds
    figure;
    plot(x_axis, young_b4cue_mean(p, :), '--', 'color', [0 0.4470 0.7410],  'LineWidth',1.5) %  (light blue)
    hold on
    jbfill(x_axis, [young_b4cue_mean(p, :) + young_b4cue_se(p, :)], [young_b4cue_mean(p, :) - young_b4cue_se(p, :)], [0 0.4470 0.7410],[0 0.4470 0.7410], 1, 0.2)
    hold on
    plot(x_axis, young_aftercue_mean(p, :), 'color', 'k',  'LineWidth',1.5)
    hold on
    jbfill(x_axis, [young_aftercue_mean(p, :) + young_aftercue_se(p, :)], [young_aftercue_mean(p, :) - young_aftercue_se(p, :)], 'k', 'k', 1, 0.2)
    hold on
    plot(x_axis, zeros(length(x_axis), 1), ':k')
    hold off
    box off
    ax = gca;
    if t == 1
        axis([-inf inf -2 2])
    else
        axis([-inf inf -2 2])
    end
    ax.LineWidth = 2.5;
    ax.FontSize = 28;
    ax.FontName = 'Arial';
    xlabel('Time (s)', 'FontSize', 32, 'FontWeight','normal')
    ylabel('Amplitude (\muV)', 'FontSize', 32, 'FontWeight','normal')
    name_channels =  channels.chanlocs(p).labels;
    title(['Young - ', name_channels] , 'FontSize', 32, 'FontWeight','normal')
%         basefilename2 = sprintf('Tepochplot_gng-sRT_young_cut_%d.fig', p);
%         savefig(basefilename2);


    figure;
    plot(x_axis, older_b4cue_mean(p, :), '--', 'color', [0 0.4470 0.7410],  'LineWidth',1.5) %plot mean gng (light blue)
    hold on
    jbfill(x_axis, [older_b4cue_mean(p, :) + older_b4cue_se(p, :)], [older_b4cue_mean(p, :) - older_b4cue_se(p, :)], [0 0.4470 0.7410],[0 0.4470 0.7410], 1, 0.2)
    hold on
    plot(x_axis, older_aftercue_mean(p, :), 'color', 'k',  'LineWidth',1.5) %plot mean gng (light blue)
    hold on
    jbfill(x_axis, [older_aftercue_mean(p, :) + older_aftercue_se(p, :)], [older_aftercue_mean(p, :) - older_aftercue_se(p, :)], 'k', 'k', 1, 0.2)
    hold on
    plot(x_axis, zeros(length(x_axis), 1), ':k')
    hold off
    box off
    ax = gca;
    if t == 1
        axis([-inf inf -2 2])
    else
        axis([-inf inf -2 2])
    end
    ax.LineWidth = 2.5;
    ax.FontSize = 28;
    ax.FontName = 'Arial';
    xlabel('Time (s)', 'FontSize', 32, 'FontWeight','normal')
    ylabel('Amplitude (\muV)', 'FontSize', 32, 'FontWeight','normal')
    name_channels =  channels.chanlocs(p).labels;
    title(['Older - ', name_channels] , 'FontSize', 32, 'FontWeight','normal')

end  


%% plot data
tasks = {'passive', 'simple RT', 'go/no-go'};
for t = 2:length(tasks)

% HEP before cue
young_b4cue_mean = squeeze(mean(HEP_T_b4cue_detrended{1, 2}(:, t, :, :), 1, 'omitnan'));
older_b4cue_mean = squeeze(mean(HEP_T_b4cue_detrended{2, 2}(:, t, :, :), 1, 'omitnan'));

young_b4cue_se = squeeze(std(HEP_T_b4cue_detrended{1, 2}(:, t, :, :), [], 1, 'omitnan'))/sqrt(size(HEP_T_b4cue_detrended{1, 2}, 1));
older_b4cue_se = squeeze(std(HEP_T_b4cue_detrended{2, 2}(:, t, :, :), [], 1, 'omitnan'))/sqrt(size(HEP_T_b4cue_detrended{1, 2}, 1));


% HEP after cue
young_aftercue_mean = squeeze(mean(HEP_T_aftercue_detrended{1, 2}(:, t, :, :), 1, 'omitnan'));
older_aftercue_mean = squeeze(mean(HEP_T_aftercue_detrended{2, 2}(:, t, :, :), 1, 'omitnan'));

young_aftercue_se = squeeze(std(HEP_T_aftercue_detrended{1, 2}(:, t, :, :), [], 1, 'omitnan'))/sqrt(size(HEP_T_aftercue_detrended{1, 2}, 1));
older_aftercue_se = squeeze(std(HEP_T_aftercue_detrended{2, 2}(:, t, :, :), [], 1, 'omitnan'))/sqrt(size(HEP_T_aftercue_detrended{1, 2}, 1));

% figure;
% topoplot(young_b4cue_mean(:, 300), channels.chanlocs); 
% Axis.FontSize = 16; caxis([-1 1]); 
% colorbar;
% % colorbar('Ticks',[0, 2, 4, 6], 'FontSize', 30, 'FontWeight','normal');
% colormap(crameri('hawaii'));
% title(['young b4cue ', tasks{t}] , 'FontSize', 32, 'FontWeight','normal')
% 
% figure;
% topoplot(young_aftercue_mean(:, 300), channels.chanlocs); 
% Axis.FontSize = 16; caxis([-1 1]); 
% colorbar;
% % colorbar('Ticks',[0, 2, 4, 6], 'FontSize', 30, 'FontWeight','normal');
% colormap(crameri('hawaii'));
% title(['young aftercue ', tasks{t}] , 'FontSize', 32, 'FontWeight','normal')

% close all
% channels = F3,Fz,F4,C3,Cz,C4,P3,Pz,P4
channels_num = 34;%[7, 16, 25, 34, 51];

    for p = channels_num

        x_axis=-199:2:400; % in mseconds
        figure;
        plot(x_axis, young_b4cue_mean(p, :), '--', 'color', [0 0.4470 0.7410],  'LineWidth',1.5) %plot mean gng (light blue)
        hold on
        jbfill(x_axis, [young_b4cue_mean(p, :) + young_b4cue_se(p, :)], [young_b4cue_mean(p, :) - young_b4cue_se(p, :)], [0 0.4470 0.7410],[0 0.4470 0.7410], 1, 0.2)
        hold on
        plot(x_axis, young_aftercue_mean(p, :), 'color', 'k',  'LineWidth',1.5) %plot mean gng (light blue)
        hold on
        jbfill(x_axis, [young_aftercue_mean(p, :) + young_aftercue_se(p, :)], [young_aftercue_mean(p, :) - young_aftercue_se(p, :)], 'k', 'k', 1, 0.2)
        hold on
        plot(x_axis, zeros(length(x_axis), 1), ':k')
        hold off
        box off
        ax = gca;
        if t == 1
            axis([-inf inf -1.5 1.5])
        else
            axis([-inf inf -2 2])
        end
        ax.LineWidth = 2.5;
        ax.FontSize = 28;
        ax.FontName = 'Arial';
        xlabel('Time (ms)', 'FontSize', 32, 'FontWeight','normal')
        ylabel('Amplitude (\muV)', 'FontSize', 32, 'FontWeight','normal')
        name_channels =  channels.chanlocs(p).labels;
        title(['young - ', tasks{t}, ' - ', name_channels] , 'FontSize', 32, 'FontWeight','normal')
%         legend({'before cue', 'after cue'}, 'Location','southeastoutside')
%         basefilename2 = sprintf('Tepochplot_gng-sRT_young_cut_%d.fig', p);
%         savefig(basefilename2);


        figure;
        plot(x_axis, older_b4cue_mean(p, :), '--', 'color', [0 0.4470 0.7410],  'LineWidth',1.5) %plot mean gng (light blue)
        hold on
        jbfill(x_axis, [older_b4cue_mean(p, :) + older_b4cue_se(p, :)], [older_b4cue_mean(p, :) - older_b4cue_se(p, :)], [0 0.4470 0.7410],[0 0.4470 0.7410], 1, 0.2)
        hold on
        plot(x_axis, older_aftercue_mean(p, :), 'color', 'k',  'LineWidth',1.5) %plot mean gng (light blue)
        hold on
        jbfill(x_axis, [older_aftercue_mean(p, :) + older_aftercue_se(p, :)], [older_aftercue_mean(p, :) - older_aftercue_se(p, :)], 'k', 'k', 1, 0.2)
        hold on
        plot(x_axis, zeros(length(x_axis), 1), ':k')
        hold off
        box off
        ax = gca;
        if t == 1
            axis([-inf inf -1.5 1.5])
        else
            axis([-inf inf -2 2])
        end
        ax.LineWidth = 2.5;
        ax.FontSize = 28;
        ax.FontName = 'Arial';
        xlabel('Time (ms)', 'FontSize', 32, 'FontWeight','normal')
        ylabel('Amplitude (\muV)', 'FontSize', 32, 'FontWeight','normal')
        name_channels =  channels.chanlocs(p).labels;
        title(['Older - ', tasks{t}, ' - ', name_channels] , 'FontSize', 32, 'FontWeight','normal')
%         legend({'before cue', 'after cue'}, 'Location','southeastoutside')


        figure;
        plot(x_axis, older_b4cue_mean(p, :), '--', 'color', 'r',  'LineWidth',1.5) %plot mean gng (light blue)
        hold on
        jbfill(x_axis, [older_b4cue_mean(p, :) + older_b4cue_se(p, :)], [older_b4cue_mean(p, :) - older_b4cue_se(p, :)], 'r', 'r', 1, 0.2)
        hold on
        plot(x_axis, young_b4cue_mean(p, :), 'color', 'k',  'LineWidth',1.5) %plot mean gng (light blue)
        hold on
        jbfill(x_axis, [young_b4cue_mean(p, :) + young_b4cue_se(p, :)], [young_b4cue_mean(p, :) - young_b4cue_se(p, :)], 'k', 'k', 1, 0.2)
        hold on
        plot(x_axis, zeros(length(x_axis), 1), ':k')
        hold off
        box off
        ax = gca;
        if t == 1
            axis([-inf inf -1.5 1.5])
        else
            axis([-inf inf -2 2])
        end
        ax.LineWidth = 2.5;
        ax.FontSize = 28;
        ax.FontName = 'Arial';
        xlabel('Time (ms)', 'FontSize', 32, 'FontWeight','normal')
        ylabel('Amplitude (\muV)', 'FontSize', 32, 'FontWeight','normal')
        name_channels =  channels.chanlocs(p).labels;
        title(['b4cue ', tasks{t}, name_channels] , 'FontSize', 32, 'FontWeight','normal')
        
        
        figure;
        plot(x_axis, older_aftercue_mean(p, :), '--', 'color', 'r',  'LineWidth',1.5) %plot mean gng (light blue)
        hold on
        jbfill(x_axis, [older_aftercue_mean(p, :) + older_aftercue_se(p, :)], [older_aftercue_mean(p, :) - older_aftercue_se(p, :)], 'r', 'r', 1, 0.2)
        hold on
        plot(x_axis, young_aftercue_mean(p, :), 'color', 'k',  'LineWidth',1.5) %plot mean gng (light blue)
        hold on
        jbfill(x_axis, [young_aftercue_mean(p, :) + young_aftercue_se(p, :)], [young_aftercue_mean(p, :) - young_aftercue_se(p, :)], 'k', 'k', 1, 0.2)
        hold on
        plot(x_axis, zeros(length(x_axis), 1), ':k')
        hold off
        box off
        ax = gca;
        if t == 1
            axis([-inf inf -1.5 1.5])
        else
            axis([-inf inf -2 2])
        end
        ax.LineWidth = 2.5;
        ax.FontSize = 28;
        ax.FontName = 'Arial';
        xlabel('Time (ms)', 'FontSize', 32, 'FontWeight','normal')
        ylabel('Amplitude (\muV)', 'FontSize', 32, 'FontWeight','normal')
        name_channels =  channels.chanlocs(p).labels;
        title(['aftercue ', tasks{t}, name_channels] , 'FontSize', 32, 'FontWeight','normal')
%         legend({'before cue', 'after cue'}, 'Location','southeastoutside')
    end
end


%% plot before cue vs after cue - both groups together - simple RT and gng
% tasks together
% HEP_T_b4cue_detrended{grp, 2}(p, t, :, :) = participants x task x
% channels x time

young_data = squeeze(mean(HEP_T_b4cue_detrended{1, 2}(:, 2:3, :, :), 2, 'omitnan'));
older_data = squeeze(mean(HEP_T_b4cue_detrended{2, 2}(:, 2:3, :, :), 2, 'omitnan'));

data = cat(1, young_data, older_data);

% HEP before cue
b4cue_mean = squeeze(mean(data, 1, 'omitnan'));
b4cue_se = squeeze(std(data, [], 1, 'omitnan'))/sqrt(size(data, 1));

% HEP after cue
young_data = squeeze(mean(HEP_T_aftercue_detrended{1, 2}(:, 2:3, :, :), 2, 'omitnan'));
older_data = squeeze(mean(HEP_T_aftercue_detrended{2, 2}(:, 2:3, :, :), 2, 'omitnan'));

data = cat(1, young_data, older_data);

% HEP before cue
aftercue_mean = squeeze(mean(data, 1, 'omitnan'));
aftercue_se = squeeze(std(data, [], 1, 'omitnan'))/sqrt(size(data, 1));


channels_num = [7, 16, 25, 34, 51];

for p = channels_num

    x_axis=-299:2:400; % in mseconds
    figure;
    plot(x_axis, b4cue_mean(p, :), 'color', [0 0.4470 0.7410],  'LineWidth',1.5) %plot mean gng (light blue)
    hold on
    jbfill(x_axis, [b4cue_mean(p, :) + b4cue_se(p, :)], [b4cue_mean(p, :) - b4cue_se(p, :)], [0 0.4470 0.7410],[0 0.4470 0.7410], 1, 0.2)
    hold on
    plot(x_axis, aftercue_mean(p, :), 'color', 'k',  'LineWidth',1.5) %plot mean gng (light blue)
    hold on
    jbfill(x_axis, [aftercue_mean(p, :) + aftercue_se(p, :)], [aftercue_mean(p, :) - aftercue_se(p, :)], 'k', 'k', 1, 0.2)
    hold on
    plot(x_axis, zeros(length(x_axis), 1), ':k')
    hold off
    box off
    ax = gca;
    if t == 1
        axis([-inf inf -1.5 1.5])
    else
        axis([-inf 400 -2.5 1.5])
    end
    ax.LineWidth = 2.5;
    ax.FontSize = 28;
    ax.FontName = 'Arial';
    xlabel('Time (ms)', 'FontSize', 32, 'FontWeight','normal')
    ylabel('Amplitude (\muV)', 'FontSize', 32, 'FontWeight','normal')
    name_channels =  channels.chanlocs(p).labels;
    title(name_channels, 'FontSize', 32, 'FontWeight','normal')
%     legend(' before cue', 'b4', 'after cue', 'after')
end


%% plot simple RT vs gng - both groups together - b4cue and after cue together
% tasks together
% HEP_T_b4cue_detrended{grp, 2}(p, t, :, :) = participants x task x
% channels x time
simpleRT_data = []; gng_data =[];
simpleRT_data(:, 1, :, :) = cat(1, HEP_T_b4cue_detrended{1, 2}(:, 2, :, :), HEP_T_b4cue_detrended{2, 2}(:, 2, :, :));
simpleRT_data(:, 2, :, :) = cat(1, HEP_T_aftercue_detrended{1, 2}(:, 2, :, :), HEP_T_aftercue_detrended{2, 2}(:, 2, :, :));
simpleRT_data = squeeze(mean(simpleRT_data, 2));

gng_data(:, 1, :, :) = cat(1, HEP_T_b4cue_detrended{1, 2}(:, 3, :, :), HEP_T_b4cue_detrended{2, 2}(:, 3, :, :));
gng_data(:, 2, :, :) = cat(1, HEP_T_aftercue_detrended{1, 2}(:, 3, :, :), HEP_T_aftercue_detrended{2, 2}(:, 3, :, :));
gng_data = squeeze(mean(gng_data, 2));

% HEP simpleRT
simpleRT_mean = squeeze(mean(simpleRT_data, 1, 'omitnan'));
simpleRT_se = squeeze(std(simpleRT_data, [], 1, 'omitnan'))/sqrt(size(simpleRT_data, 1));

% HEP gng
gng_mean = squeeze(mean(gng_data, 1, 'omitnan'));
gng_se = squeeze(std(gng_data, [], 1, 'omitnan'))/sqrt(size(gng_data, 1));

channels_num = [7, 16, 25, 34, 51];

for p = channels_num

    x_axis=-199:2:400; % in mseconds
    figure;
    plot(x_axis, simpleRT_mean(p, :), 'color', [0 0.4470 0.7410],  'LineWidth',1.5) %plot mean gng (light blue)
    hold on
    jbfill(x_axis, [simpleRT_mean(p, :) + simpleRT_se(p, :)], [simpleRT_mean(p, :) - simpleRT_se(p, :)], [0 0.4470 0.7410],[0 0.4470 0.7410], 1, 0.2)
    hold on
    plot(x_axis, gng_mean(p, :), 'color', 'k',  'LineWidth',1.5) %plot mean gng (light blue)
    hold on
    jbfill(x_axis, [gng_mean(p, :) + gng_se(p, :)], [gng_mean(p, :) - gng_se(p, :)], 'k', 'k', 1, 0.2)
    hold on
    plot(x_axis, zeros(length(x_axis), 1), ':k')
    hold off
    box off
    ax = gca;
    if t == 1
        axis([-inf inf -1.5 1.5])
    else
        axis([-inf inf -2.5 1.5])
    end
    ax.LineWidth = 2.5;
    ax.FontSize = 28;
    ax.FontName = 'Arial';
    xlabel('Time (ms)', 'FontSize', 32, 'FontWeight','normal')
    ylabel('Amplitude (\muV)', 'FontSize', 32, 'FontWeight','normal')
    name_channels =  channels.chanlocs(p).labels;
    title(name_channels, 'FontSize', 32, 'FontWeight','normal')
%     legend(' before cue', 'b4', 'after cue', 'after')
end

