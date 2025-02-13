% plot HEP - all time course
% all participants together and group differences
% exclude participants that were not included in T peak analyses: 21 50 55
% - T peak not detected
clear; close all
load subjects
load HEP4plot_avg; % calculated in \LIMO_stats\HEP_grp_task_effect_RLocked_nobaseline_filt\LIMO_HEP_grp_task_effect_Rlocked_nobaseline.m
%  HEP4plot_avg{group, t}(yng, :, :) = mean(EEG.data, 3);
% delete excluded participants
older = [11 12 14 20 21 22 23 32 37 38 41 43 47 48 49 52 55 57 58 63 64 65 67 69 7 70 71 75 8 83 86];
young = [13 15 16 25 26 28 31 33 34 36 4 42 44 45 46 50 51 53 54 56 59 6 62 66 68 72 74 76 78 80 82 84 85 9];



subj_older = []; subj_young = [];
for s = 1:length(subjects)
    if strcmp(subjects{s}(4), '_')
        sbj = 3;
    else
        sbj = 3:4;
    end
    subj_number = str2num(subjects{s}(sbj));
    if ismember(subj_number, young)
        subj_young = [subj_young; subj_number];
    else
        subj_older = [subj_older; subj_number];
    end
end

exc_young = find(subj_young == 50)
exc_older = ismember(subj_older, [21, 55])

HEP4plot_avg{1, 1}(exc_young, :, :) = [];
HEP4plot_avg{1, 2}(exc_young, :, :) = [];

HEP4plot_avg{2, 1}(exc_older, :, :) = [];
HEP4plot_avg{2, 2}(exc_older, :, :) = [];

load('E:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\LIMO_stats\expected_chanlocs_both.mat');

% average across tasks
for grp = 1:2
    for prt = 1:size(HEP4plot_avg{grp, 1}, 1)
        clear tmp
        tmp(:, :, 1) = squeeze(HEP4plot_avg{grp, 1}(prt, :, :));
        tmp(:, :, 2) = squeeze(HEP4plot_avg{grp, 2}(prt, :, :));
        HEP{grp}(prt, :, :) = mean(tmp, 3);
    end
end

%% average across groups - effect of task
HEP_task = cell(2, 1);

HEP_task{1} = cat(1, HEP4plot_avg{1, 1}, HEP4plot_avg{2, 1});
HEP_task{2} = cat(1, HEP4plot_avg{1, 2}, HEP4plot_avg{2, 2});

simpleRT_mean = squeeze(mean(HEP_task{1}, 1));
simpleRT_se = squeeze(std(HEP_task{1}, [], 1))/sqrt(size(HEP_task{1}, 1));

gng_mean = squeeze(mean(HEP_task{2}, 1));
gng_se = squeeze(std(HEP_task{2}, [], 1))/sqrt(size(HEP_task{2}, 1));

% color for each task
clr = [230, 159, 0; 0 114 178; 213, 94, 0]./255;
x_axis = -149:2:700; % in ms
channels_num = [16 34 51];%[7, 25, 43, 57];
for chn = channels_num

   

    figure;
    
    % grey background in time window not used for stats
    for idx = -145:5:3001
        plot([idx idx], [-1.05 1], 'color', [.9 .9 .9], 'LineWidth', 5)
        hold on
    end
    
    plot(x_axis, simpleRT_mean(chn, :), 'color', clr(2, :),  'LineWidth',1.5) %plot(light blue)
    hold on
    jbfill(x_axis, [simpleRT_mean(chn, :) + simpleRT_se(chn, :)], [simpleRT_mean(chn, :) - simpleRT_se(chn, :)],clr(2, :),clr(2, :), 1, 0.2)
    hold on
%     plot(x_axis(1:125), simpleRT_mean(chn, 1:125), 'color', clr(2, :)+.15,  'LineWidth',1.5) %plot(light blue)
%     hold on
    plot(x_axis, gng_mean(chn, :), '--', 'color', clr(3, :),  'LineWidth',1.5) 
    hold on
    jbfill(x_axis, [gng_mean(chn, :) + gng_se(chn, :)], [gng_mean(chn, :) - gng_se(chn, :)], clr(3, :), clr(3, :), 1, 0.2)
    hold on
%     plot(x_axis(1:125), gng_mean(chn, 1:125), '--', 'color', clr(3, :)+.15, 'LineWidth',1.5) 
%     hold on
    plot(x_axis, zeros(length(x_axis), 1), '--k',  'LineWidth', 1) 
    
    
    hold off
    box off
    ax = gca;
%     axis([-inf 350 -1.1 1])
    ax.LineWidth = 2.5;
    ax.FontSize = 28;
    ax.FontName = 'Arial';
    xlabel('Time (ms)', 'FontSize', 32, 'FontWeight','normal')
    ylabel('Amplitude (\muV)', 'FontSize', 32, 'FontWeight','normal')
    name_channels =  expected_chanlocs(chn).labels;
    title(name_channels , 'FontSize', 32, 'FontWeight','normal')
        
end


%% plot HEP both groups together, avg both tasks

HEP_all = cat(1, HEP{1}, HEP{2});
HEP_all_mean = squeeze(mean(HEP_all, 1));
HEP_all_se = squeeze(std(HEP_all, [], 1)/sqrt(size(HEP_all, 1)));
x_axis = -149:2:700; % in ms
level2_onesample_dir = [pwd, '\level2_onesample_HEP_grps_together'];
channels_num = [16 34 51];%[7, 25, 43, 57];
for chn = channels_num
    figure;
    plot(x_axis, HEP_all_mean(chn, :), 'color', 'k',  'LineWidth',1.5)
    hold on
    jbfill(x_axis, [HEP_all_mean(chn, :) + HEP_all_se(chn, :)], [HEP_all_mean(chn, :) - HEP_all_se(chn, :)], 'k', 'k', 1, 0.2)
    hold on
    plot(x_axis(1:125), HEP_all_mean(chn, 1:125), 'color', [.5 .5 .5],  'LineWidth',1.5)
    hold on
    plot(x_axis, zeros(length(x_axis), 1), '--k', 'LineWidth', 1)
    hold on

    %     % horizontal line marking where young group is different from zero
%     %     % find mask for channel chosen   
%     load([level2_onesample_dir, filesep, 'mask.mat'])
%     hold on
%     % plot horizontal line where the groups differ - mask starts on .1
%     % after R peak = .250 ms = 125 points
%     plot(x_axis(mask(chn, :) > 0) + 250, ones(length(find(mask(chn, :) > 0)), 1) * -1.15, 's', 'MarkerSize', 4.5, 'color', 'k', ...
%         'MarkerFaceColor', 'k')

    hold off
    box off
    ax = gca;
    axis([-inf inf -1.5 2.25])
    ax.LineWidth = 2.5;
    ax.FontSize = 28;
    ax.FontName = 'Arial';
    xlabel('Time (ms)', 'FontSize', 32, 'FontWeight','normal')
    ylabel('Amplitude (\muV)', 'FontSize', 32, 'FontWeight','normal')
    name_channels =  expected_chanlocs(chn).labels;
    title(name_channels , 'FontSize', 32, 'FontWeight','normal')
end

%% calculate average amplitude inside clusters where HEP is significantly different from zero
% load([pwd, '\level2_onesample_HEP_grps_together\mask.mat'])
% 
% [row, column] = find(mask == 1);
% 
% avg_ampl_cluster1 = mean(mean(HEP_all_mean(row, column + 125), 1));
% 
% % across participants standard deviation
% HEP_all_std = squeeze(std(HEP_all, [], 1));
% std_ampl_cluster1 = mean(mean(HEP_all_std(row, column + 125), 1));

%% calculate average amplitude inside clusters - older group - one-sample t-test analysis

% [row, column] = find(mask == 1);
% older_all_mean = squeeze(mean(HEP{2}(:, 1:59, :), 1));
% avg_ampl_cluster1 = mean(mean(older_all_mean(row, column + 125), 1))

%% calculate average group difference in cluster presenting an effect of group
% load('E:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\LIMO_stats\HEP_grp_task_effect_Tlocked\level2_RM_ANOVA_HEP_grp_task_effect_50-400ms\mask_group_effect')
%  
% [row, column] = find(mask == 1);
% older_all_mean = squeeze(mean(HEP{2}(:, 1:59, :), 1));
% avg_ampl_cluster1 = mean(mean(older_all_mean(row, column + 125), 1))
% 
% young_all_mean = squeeze(mean(HEP{1}(:, 1:59, :), 1));
% avg_ampl_cluster2 = mean(mean(young_all_mean(row, column + 125), 1))
% 
% avg_ampl_cluster1 - avg_ampl_cluster2


%%
young_all_mean = squeeze(mean(HEP{1}, 1));
young_all_se = squeeze(std(HEP{1}, [], 1))/sqrt(size(HEP{1}, 1));

older_all_mean = squeeze(mean(HEP{2}, 1));
older_all_se = squeeze(std(HEP{2}, [], 1))/sqrt(size(HEP{2}, 1));

% comparison across groups - including time windows with group differences
% level2_RM_ANOVA_dir = [pwd filesep 'level2_RM_ANOVA_HEP_grp_task_effect_50-400ms'];
% load([level2_RM_ANOVA_dir, filesep, 'mask_group_effect.mat'])

% calculate amplitude difference across groups in cluster showing group difference

% [row, column] = find(mask == 1);
% 
% avg_ampl_cluster_young = mean(mean(young_all_mean(row, column + 125), 1));
% avg_ampl_cluster_older = mean(mean(older_all_mean(row, column + 125), 1));
% 
% avg_ampl_cluster = avg_ampl_cluster_older-avg_ampl_cluster_young


x_axis = -149:2:700; % in ms
channels_num = [16 34 51];%[7, 25, 43, 57];
for chn = channels_num

    figure;
    
    % grey background in time window not used for stats
    for idx = -135:5:325
        plot([idx idx], [-1.15 2], 'color', [.9 .9 .9], 'LineWidth', 5)
        hold on
    end
    plot(x_axis, young_all_mean(chn, :), 'color', 'k',  'LineWidth',1.5)
    hold on
    jbfill(x_axis, [young_all_mean(chn, :) + young_all_se(chn, :)], [young_all_mean(chn, :) - young_all_se(chn, :)], 'k', 'k', 1, 0.2)
    hold on
%     plot(x_axis(1:125), young_all_mean(chn, 1:125), 'color', [.6 .6 .6],  'LineWidth',1.5) %plot(light blue)
%     hold on
    plot(x_axis, older_all_mean(chn, :), '--r',  'LineWidth',1.5) 
    hold on
    jbfill(x_axis, [older_all_mean(chn, :) + older_all_se(chn, :)], [older_all_mean(chn, :) - older_all_se(chn, :)], 'r', 'r', 1, 0.2)
    hold on
%     plot(x_axis(1:125), older_all_mean(chn, 1:125), '--', 'color', [1 .6 .6], 'LineWidth',1.5) 
%     hold on
    plot(x_axis, zeros(length(x_axis), 1), '--k',  'LineWidth', 1.5) 
    hold on
    plot([0 0], [-1.2 2], '--k',  'LineWidth', 1.5) 
    hold on
    
%     % horizontal line marking where young group is different from zero
%     % find mask for channel chosen   
%     level2_onesample_dir = [pwd, '\level2_onesample_HEP_separate_grps_330-630'];
%     if exist([level2_onesample_dir, filesep, 'young', filesep, 'mask.mat'], 'file')
%         load([level2_onesample_dir, filesep, 'young', filesep, 'mask.mat'])
%         plot(x_axis(find(mask(chn, :) > 0) + 240), ones(length(find(mask(chn, :) > 0)), 1) * -.4, 'color', [.5 .5 .5], 'LineWidth', 4.5)
%         hold on
%     end
%     if exist([level2_onesample_dir, filesep, 'older', filesep, 'mask.mat'], 'file')
%         load([level2_onesample_dir, filesep, 'older', filesep, 'mask.mat'])
%         plot(x_axis(find(mask(chn, :) > 0) + 240), ones(length(find(mask(chn, :) > 0)), 1) * -.6, 's', 'MarkerSize', 4.5, 'color', [1 .5 .5], ...
%             'MarkerFaceColor', [1 .5 .5])
%         hold on
%     end
    
    % horizontal line marking where the HEP is different from zero - all
    % groups together
    % find mask for channel chosen 
    level2_onesample_dir = [pwd, filesep, 'level2_onesample_HEP_grpstogeth_330-630ms_Tsubjs'];
    if exist([level2_onesample_dir, filesep, 'mask.mat'], 'file')
        load([level2_onesample_dir, filesep, 'mask.mat'])
        plot(x_axis(find(mask(chn, :) > 0) + 240), ones(length(find(mask(chn, :) > 0)), 1) * -.6, 's', 'MarkerSize', 5.5, 'color', 'm', ...
            'MarkerFaceColor', 'm')
        hold on
    end
    
    
    % plot horizontal line where the groups differ - mask starts on .1
    % after R peak = .3 s = 225 points
    level2_RM_ANOVA_dir = [pwd, '\level2_RM_ANOVA_HEP_grp_task_effect_330-630_Tsubjs'];
    if exist([level2_RM_ANOVA_dir, filesep 'mask_group_effect.mat'], 'file')
        load([level2_RM_ANOVA_dir, filesep 'mask_group_effect.mat'])
        plot(x_axis(find(mask(chn, :) > 0) + 240) , ones(length(find(mask(chn, :) > 0)), 1) * -.9, 's', 'MarkerSize', 5.5, 'color', 'k', ...
        'MarkerFaceColor', 'k')
    end
    
    hold off
    box off
    ax = gca;
    axis([-inf 630 -1.2 2])
    ax.LineWidth = 2.5;
    ax.FontSize = 28;
    ax.FontName = 'Arial';
    xlabel('Time (ms)', 'FontSize', 32, 'FontWeight','normal')
    ylabel('Amplitude (\muV)', 'FontSize', 32, 'FontWeight','normal')
    name_channels =  expected_chanlocs(chn).labels;
    title(name_channels , 'FontSize', 32, 'FontWeight','normal')
        
end

%% HEP amplitude difference across groups in significant mask
% calculate average amplitude difference across groups inside clusters where there is a significant task effect
level2_RM_ANOVA_dir = [pwd, '\level2_RM_ANOVA_HEP_grp_task_effect_330-630_Tsubjs'];
load([level2_RM_ANOVA_dir, filesep, 'mask_group_effect.mat'])

% mask starts at 330 ms - HEP starts at -149 = 330+149 = 479/2 = 240 data points

[row, column] = find(mask > 0);
young_mask_mean = [];
for s = 1:size(HEP{1}, 1)
    tmp = [];
    for idx = 1:length(row)
        tmp = [tmp; HEP{1}(s, row(idx), column(idx) + 240)];
    end
    young_mask_mean(s) = mean(tmp);
end

young_mean = mean(young_mask_mean)
young_std = std(young_mask_mean)

older_mask_mean = [];
for s = 1:size(HEP{2}, 1)
    tmp = [];
    for idx = 1:length(row)
        tmp = [tmp; HEP{2}(s, row(idx), column(idx) + 240)];
    end
    older_mask_mean(s) = mean(tmp);
end
        
older_mean = mean(older_mask_mean)
older_std = std(older_mask_mean)




% std_ampl_cluster_young = mean(mean(young_std(row, column + 240), 1)) %.6202
% std_ampl_cluster_older = mean(mean(older_std(row, column + 240), 1)) % .5983


std_ampl_cluster_young = std(mean(mean(HEP{1}(:, row, column + 240), 3), 2)) %.6202
std_ampl_cluster_older = std(mean(mean(HEP{2}(:, row, column + 240), 3), 2)) % .5983


avg_ampl_cluster_older-avg_ampl_cluster_young

HEP_groupdiff_mean = squeeze(mean(HEP{2}, 1)) - squeeze(mean(HEP{1}, 1));

groupdiff_mean = mean(mean(HEP_groupdiff_mean(row, column + 240), 1))




%% data for topoplot - difference across groups
% find channels where task effect is significant
level2_RM_ANOVA_dir = [pwd, '\level2_RM_ANOVA_HEP_grp_task_effect_330-630_Tsubjs'];
load([level2_RM_ANOVA_dir, filesep, 'mask_group_effect.mat'])
% average across time window where task effect is significant
time_window = find(sum(mask, 1) > 0);

young_topoplot = mean(young_all_mean(:, 240+time_window), 2);
older_topoplot = mean(older_all_mean(:, 240+time_window), 2);


% channels with significant group difference
sig_chan_number = find(sum(mask, 2) > 0);


figure;
topoplot(ones(59, 1)*.025, expected_chanlocs, 'electrodes', 'off', 'plotchans', sig_chan_number, 'style', 'blank',...
    'plotdisk', 'on',  'hcolor'  , 'none') ; hold on
topoplot(young_topoplot, expected_chanlocs, 'electrodes', 'off'); 
% caxis([-.5 .5]); %c.Axis.FontSize = 16;
cH = colorbar; set(cH,'FontSize',30);
% colorbar('Ticks',[0, 1, 2, 3, 4, 5], 'FontSize', 18, 'FontWeight','normal');
colormap(crameri('vik')); % needs colour maps from http://www.fabiocrameri.ch/colourmaps.php
% title('Young', 'FontSize', 30, 'FontWeight','normal')
% set(get(gca,'title'),'Position',[0,-.65, 0])
%         text(5, 0.4, title_txt{task, grp})


figure;
topoplot(ones(59, 1)*.025, expected_chanlocs, 'electrodes', 'off', 'plotchans', sig_chan_number, 'style', 'blank',...
    'plotdisk', 'on',  'hcolor'  , 'none') ; hold on
topoplot(older_topoplot, expected_chanlocs, 'electrodes', 'off'); 
% caxis([-.5 .5]); %c.Axis.FontSize = 16;
cH = colorbar; set(cH,'FontSize',30);
% colorbar('Ticks',[0, 1, 2, 3, 4, 5], 'FontSize', 18, 'FontWeight','normal');
colormap(crameri('vik')); % needs colour maps from http://www.fabiocrameri.ch/colourmaps.php
% title('Older', 'FontSize', 30, 'FontWeight','normal')
% set(get(gca,'title'),'Position',[0,-.65, 0])


figure;
topoplot(ones(59, 1)*.025, expected_chanlocs, 'electrodes', 'off', 'plotchans', sig_chan_number, 'style', 'blank',...
    'plotdisk', 'on',  'hcolor'  , 'none') ; hold on
topoplot(older_topoplot-young_topoplot, expected_chanlocs, 'electrodes', 'off'); 
% caxis([-.4 .4]); %c.Axis.FontSize = 16;
cH = colorbar; set(cH,'FontSize',30);
% colorbar('Ticks',[0, 1, 2, 3, 4, 5], 'FontSize', 18, 'FontWeight','normal');
colormap(crameri('vik')); % needs colour maps from http://www.fabiocrameri.ch/colourmaps.php



%% plot data using trimmed means
% [t,tmdata,trimci,se,p,tcrit,df]=limo_trimci(data, percent, alphav, nullvalue)

for s=1:size(HEP{1}, 1)
    data_young(:,:, s) = HEP{1}(s,:, :);
end

for s=1:size(HEP{2}, 1)
    data_older(:,:, s) = HEP{2}(s,:, :);
end

[~,young_all_mean,~,young_all_se,~,~,~]=limo_trimci(data_young); % 20% trimmed means

[~,older_all_mean,~,older_all_se,~,~,~]=limo_trimci(data_older); % 20% trimmed means


channels_num = [16 34 51];%[7, 25, 43, 57];
for chn = channels_num

   
    x_axis = -199:2:400; % in ms
    figure;
    
    % grey background in time window not used for stats
    for idx = -185:5:50
        plot([idx idx], [-1.15 1], 'color', [.9 .9 .9], 'LineWidth', 5)
        hold on
    end
    
    plot(x_axis, young_all_mean(chn, :), 'color', 'k',  'LineWidth',1.5)
    hold on
    jbfill(x_axis, [young_all_mean(chn, :) + young_all_se(chn, :)], [young_all_mean(chn, :) - young_all_se(chn, :)], 'k', 'k', 1, 0.2)
    hold on
%     plot(x_axis(1:125), young_all_mean(chn, 1:125), 'color', [.6 .6 .6],  'LineWidth',1.5) %plot(light blue)
%     hold on
    plot(x_axis, older_all_mean(chn, :), '--r',  'LineWidth',1.5) 
    hold on
    jbfill(x_axis, [older_all_mean(chn, :) + older_all_se(chn, :)], [older_all_mean(chn, :) - older_all_se(chn, :)], 'r', 'r', 1, 0.2)
    hold on
%     plot(x_axis(1:125), older_all_mean(chn, 1:125), '--', 'color', [1 .6 .6], 'LineWidth',1.5) 
%     hold on
    plot(x_axis, zeros(length(x_axis), 1), '--k',  'LineWidth', 1) 

    hold off
    box off
    ax = gca;
    axis([-inf 350 -1.2 1])
    ax.LineWidth = 2.5;
    ax.FontSize = 28;
    ax.FontName = 'Arial';
    xlabel('Time (ms)', 'FontSize', 32, 'FontWeight','normal')
    ylabel('Amplitude (\muV)', 'FontSize', 32, 'FontWeight','normal')
    name_channels =  expected_chanlocs(chn).labels;
    title(name_channels , 'FontSize', 32, 'FontWeight','normal')
        
end


%% data for topoplot
% find channels where HEP is significantly different from zero
level2_RM_ANOVA_dir = [pwd filesep 'level2_onesample_HEP_grps_together'];
load([level2_RM_ANOVA_dir, filesep, 'mask.mat'])
% average across time window where HEP is significantly different from zero
time_window = find(sum(mask, 1) > 0);

young_topoplot = mean(young_all_mean(:, 125+time_window), 2);
older_topoplot = mean(older_all_mean(:, 125+time_window), 2);


% channels with significant group difference
sig_chan_number = find(sum(mask, 2) > 0);

title_txt = {'Young' 'Older'};
       
%         for x = 1:length(sig_chan_number)
%             sig_chans{grp, task}{x} = chanlocs_EEGChanOnly(sig_chan_number(x)).labels;
%         end

figure;
topoplot(ones(59, 1)*.025, expected_chanlocs, 'electrodes', 'off', 'plotchans', sig_chan_number, 'style', 'blank',...
    'plotdisk', 'on',  'hcolor'  , 'none') ; hold on
topoplot(young_topoplot, expected_chanlocs, 'electrodes', 'off'); 
caxis([-.4 .4]); %c.Axis.FontSize = 16;
cH = colorbar; set(cH,'FontSize',30);
% colorbar('Ticks',[0, 1, 2, 3, 4, 5], 'FontSize', 18, 'FontWeight','normal');
colormap(crameri('vik')); % needs colour maps from http://www.fabiocrameri.ch/colourmaps.php
% title('Young', 'FontSize', 30, 'FontWeight','normal')
% set(get(gca,'title'),'Position',[0,-.65, 0])
%         text(5, 0.4, title_txt{task, grp})


figure;
topoplot(ones(59, 1)*.025, expected_chanlocs, 'electrodes', 'off', 'plotchans', sig_chan_number, 'style', 'blank',...
    'plotdisk', 'on',  'hcolor'  , 'none') ; hold on
topoplot(older_topoplot, expected_chanlocs, 'electrodes', 'off'); 
caxis([-.4 .4]); %c.Axis.FontSize = 16;
cH = colorbar; set(cH,'FontSize',30);
% colorbar('Ticks',[0, 1, 2, 3, 4, 5], 'FontSize', 18, 'FontWeight','normal');
colormap(crameri('vik')); % needs colour maps from http://www.fabiocrameri.ch/colourmaps.php
% title('Older', 'FontSize', 30, 'FontWeight','normal')
% set(get(gca,'title'),'Position',[0,-.65, 0])


%% plot topoplots - all groups all tasks together

mean_data = squeeze(mean(HEP_all(:, 1:59, 175:225), 1));
all_topoplot = mean(mean_data, 2);

% channels where the T-HEP is significantly different from zero
load([pwd, '/level2_onesample_HEP_all_partcts_50-400ms/mask']) 
sig_chan_number = find(sum(mask, 2) > 0);

title_txt = {'Young' 'Older'};
       
%         for x = 1:length(sig_chan_number)
%             sig_chans{grp, task}{x} = chanlocs_EEGChanOnly(sig_chan_number(x)).labels;
%         end

figure;
topoplot(ones(59, 1)*.025, expected_chanlocs, 'electrodes', 'off', 'plotchans', sig_chan_number, 'style', 'blank',...
    'plotdisk', 'on',  'hcolor'  , 'none') ; hold on
topoplot(all_topoplot, expected_chanlocs, 'electrodes', 'off'); 
% caxis([-1 1]); %c.Axis.FontSize = 16;
cH = colorbar; set(cH,'FontSize',30);
% colorbar('Ticks',[0, 1, 2, 3, 4, 5], 'FontSize', 18, 'FontWeight','normal');
colormap(crameri('vik')); % needs colour maps from http://www.fabiocrameri.ch/colourmaps.php
% title('Young', 'FontSize', 30, 'FontWeight','normal')
% set(get(gca,'title'),'Position',[0,-.65, 0])
%         text(5, 0.4, title_txt{task, grp})



