% plot HEP - all time course
% all participants together and group differences
clear

load HEP_T_locked; % calculated in \HEP_grp_task_effect_Tlocked\LIMO_HEPTlocked_grp_task_effect.m
% HEP_T_locked{group, t} = group x task
load('E:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\LIMO_stats\expected_chanlocs_both.mat');

% average across tasks
for grp = 1:2
    for prt = 1:size(HEP_T_locked{grp, 1}, 1)
        clear tmp
        tmp(:, :, 1) = squeeze(HEP_T_locked{grp, 1}(prt, :, :));
        tmp(:, :, 2) = squeeze(HEP_T_locked{grp, 2}(prt, :, :));
        HEP{grp}(prt, :, :) = mean(tmp, 3); 
    end
end

%% average across groups - effect of task
HEP_task = cell(2, 1);

HEP_task{1} = cat(1, HEP_T_locked{1, 1}, HEP_T_locked{2, 1});
HEP_task{2} = cat(1, HEP_T_locked{1, 2}, HEP_T_locked{2, 2});

simpleRT_mean = squeeze(mean(HEP_task{1}, 1));
simpleRT_se = squeeze(std(HEP_task{1}, [], 1))/sqrt(size(HEP_task{1}, 1));

gng_mean = squeeze(mean(HEP_task{2}, 1));
gng_se = squeeze(std(HEP_task{2}, [], 1))/sqrt(size(HEP_task{2}, 1));

% color for each task
clr = [230, 159, 0; 0 114 178; 213, 94, 0]./255;

channels_num = [16 34 51];%[7, 25, 43, 57];
for chn = channels_num

   
    x_axis = -199:2:400; % in ms
    figure;
    plot(x_axis, simpleRT_mean(chn, :), 'color', clr(2, :),  'LineWidth',1.5) %plot(light blue)
    hold on
    jbfill(x_axis, [simpleRT_mean(chn, :) + simpleRT_se(chn, :)], [simpleRT_mean(chn, :) - simpleRT_se(chn, :)],clr(2, :),clr(2, :), 1, 0.2)
    hold on
    plot(x_axis(1:125), simpleRT_mean(chn, 1:125), 'color', clr(2, :)+.15,  'LineWidth',1.5) %plot(light blue)
    hold on
    plot(x_axis, gng_mean(chn, :), '--', 'color', clr(3, :),  'LineWidth',1.5) 
    hold on
    jbfill(x_axis, [gng_mean(chn, :) + gng_se(chn, :)], [gng_mean(chn, :) - gng_se(chn, :)], clr(3, :), clr(3, :), 1, 0.2)
    hold on
    plot(x_axis(1:125), gng_mean(chn, 1:125), '--', 'color', clr(3, :)+.15, 'LineWidth',1.5) 
    hold on
    plot(x_axis, zeros(length(x_axis), 1), '--k',  'LineWidth', 1) 
    
    
    hold off
    box off
    ax = gca;
    axis([-inf inf -1.5 1.3])
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

level2_onesample_dir = [pwd, '\level2_onesample_HEP_all_partcts_50-400ms'];
channels_num = [16 34 51];%[7, 25, 43, 57];
for chn = channels_num
    x_axis = -199:2:600; % in ms
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
    %     % find mask for channel chosen   
    load([level2_onesample_dir, filesep, 'mask.mat'])
    hold on
    % plot horizontal line where the groups differ - mask starts on .1
    % after R peak = .250 ms = 125 points
    plot(x_axis(mask(chn, :) > 0) + 250, ones(length(find(mask(chn, :) > 0)), 1) * -1.15, 's', 'MarkerSize', 4.5, 'color', 'k', ...
        'MarkerFaceColor', 'k')

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

%% calculate average amplitude inside clusters

[row, column] = find(mask == 1);

avg_ampl_cluster1 = mean(mean(HEP_all_mean(row, column + 125), 1));

[row, column] = find(mask == 2);

avg_ampl_cluster2 = mean(mean(HEP_all_mean(row, column + 125), 1));

%% calculate average amplitude inside clusters - older group - one-sample t-test analysis

[row, column] = find(mask == 1);
older_all_mean = squeeze(mean(HEP{2}(:, 1:59, :), 1));
avg_ampl_cluster1 = mean(mean(older_all_mean(row, column + 125), 1))

%% calculate average group difference in cluster presenting an effect of group
load('E:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\LIMO_stats\HEP_grp_task_effect_Tlocked\level2_RM_ANOVA_HEP_grp_task_effect_50-400ms\mask_group_effect')
 
[row, column] = find(mask == 1);
older_all_mean = squeeze(mean(HEP{2}(:, 1:59, :), 1));
avg_ampl_cluster1 = mean(mean(older_all_mean(row, column + 125), 1))

young_all_mean = squeeze(mean(HEP{1}(:, 1:59, :), 1));
avg_ampl_cluster2 = mean(mean(young_all_mean(row, column + 125), 1))

avg_ampl_cluster1 - avg_ampl_cluster2


%%
young_all_mean = squeeze(mean(HEP{1}, 1));
young_all_se = squeeze(std(HEP{1}, [], 1))/sqrt(size(HEP{1}, 1));

older_all_mean = squeeze(mean(HEP{2}, 1));
older_all_se = squeeze(std(HEP{2}, [], 1))/sqrt(size(HEP{2}, 1));

% comparison across groups - including time windows with group differences
level2_RM_ANOVA_dir = [pwd filesep 'level2_RM_ANOVA_HEP_grp_task_effect_50-400ms'];
load([level2_RM_ANOVA_dir, filesep, 'mask_group_effect.mat'])

% calculate amplitude difference across groups in cluster showing group difference

[row, column] = find(mask == 1);

avg_ampl_cluster_young = mean(mean(young_all_mean(row, column + 125), 1));
avg_ampl_cluster_older = mean(mean(older_all_mean(row, column + 125), 1));

avg_ampl_cluster = avg_ampl_cluster_older-avg_ampl_cluster_young



channels_num = [16 34 51];%[7, 25, 43, 57];
for chn = channels_num

   
    x_axis = -199:2:400; % in ms
    figure;
    plot(x_axis, young_all_mean(chn, :), 'color', 'k',  'LineWidth',1.5)
    hold on
    jbfill(x_axis, [young_all_mean(chn, :) + young_all_se(chn, :)], [young_all_mean(chn, :) - young_all_se(chn, :)], 'k', 'k', 1, 0.2)
    hold on
    plot(x_axis(1:125), young_all_mean(chn, 1:125), 'color', [.6 .6 .6],  'LineWidth',1.5) %plot(light blue)
    hold on
    plot(x_axis, older_all_mean(chn, :), '--r',  'LineWidth',1.5) 
    hold on
    jbfill(x_axis, [older_all_mean(chn, :) + older_all_se(chn, :)], [older_all_mean(chn, :) - older_all_se(chn, :)], 'r', 'r', 1, 0.2)
    hold on
    plot(x_axis(1:125), older_all_mean(chn, 1:125), '--', 'color', [1 .6 .6], 'LineWidth',1.5) 
    hold on
    plot(x_axis, zeros(length(x_axis), 1), '--k',  'LineWidth', 1) 
    hold on
    
%     % horizontal line marking where young group is different from zero
%     % find mask for channel chosen   
%     level2_onesample_dir = [pwd, '\level2_onesample_HEP'];
%     load([level2_onesample_dir, filesep, 'young', filesep, 'mask.mat'])
%     plot(x_axis(mask(chn, :) == 1), ones(length(find(mask(chn, :) == 1)), 1) * -1.25, 'k', 'LineWidth', 4.5)
%     hold on
%     load([level2_onesample_dir, filesep, 'older', filesep, 'mask.mat'])
%     plot(x_axis(mask(chn, :) > 0), ones(length(find(mask(chn, :) > 0)), 1) * -1.35, 's', 'MarkerSize', 4.5, 'color', 'r', ...
%         'MarkerFaceColor', 'r')
    hold on
    % plot horizontal line where the groups differ - mask starts on .1
    % after R peak = .250 ms = 125 points
    plot(x_axis(mask(chn, :) > 0) + 250, ones(length(find(mask(chn, :) > 0)), 1) * -1.15, 's', 'MarkerSize', 4.5, 'color', 'k', ...
        'MarkerFaceColor', 'k')
    
    hold off
    box off
    ax = gca;
    axis([-inf inf -1.5 1.3])
    ax.LineWidth = 2.5;
    ax.FontSize = 28;
    ax.FontName = 'Arial';
    xlabel('Time (ms)', 'FontSize', 32, 'FontWeight','normal')
    ylabel('Amplitude (\muV)', 'FontSize', 32, 'FontWeight','normal')
    name_channels =  expected_chanlocs(chn).labels;
    title(name_channels , 'FontSize', 32, 'FontWeight','normal')
        
end

%% plot topoplots
% data for topoplot = average across time window including group difference
% = 50-150 ms after T peak

young_topoplot = mean(young_all_mean(:, 125:175), 2);
older_topoplot = mean(older_all_mean(:, 125:175), 2);

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
caxis([-.8 .8]); %c.Axis.FontSize = 16;
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
caxis([-.8 .8]); %c.Axis.FontSize = 16;
cH = colorbar; set(cH,'FontSize',30);
% colorbar('Ticks',[0, 1, 2, 3, 4, 5], 'FontSize', 18, 'FontWeight','normal');
colormap(crameri('vik')); % needs colour maps from http://www.fabiocrameri.ch/colourmaps.php
% title('Older', 'FontSize', 30, 'FontWeight','normal')
% set(get(gca,'title'),'Position',[0,-.65, 0])


%% plot topoplots - all groups all tasks together
% data for topoplot = average across time window including group difference
% = 50-150 ms after T peak

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



