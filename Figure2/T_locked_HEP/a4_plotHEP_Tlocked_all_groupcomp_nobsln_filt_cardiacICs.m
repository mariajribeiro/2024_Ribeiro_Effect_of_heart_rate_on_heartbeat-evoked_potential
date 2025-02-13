% plot HEP - all time course
% all participants together and group differences
clear; close all

load HEP_T_locked; % calculated in \HEP_grp_task_effect_Tlocked\LIMO_HEPTlocked_grp_task_effect_nobaseline.m
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

channels_num = [15 16 23 24 33 34 51];%[7, 25, 43, 57];
for chn = channels_num

   
    x_axis = -199:2:400; % in ms
    figure;
    
    % grey background in time window not used for stats
    for idx = -185:5:50
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
    hold on
    plot([0 0] , [-1.1 .7], '--k',  'LineWidth', 1) 
    hold on
    % horizontal line marking where there is an effect of task
    % find mask for channel chosen   
    level2_RM_ANOVA_dir = [pwd, '\level2_RM_ANOVA_HEP_grp_task_effect_50-350ms'];
    load([level2_RM_ANOVA_dir , filesep, 'mask_main_effect.mat'])
    plot(x_axis(find((mask(chn, :) > 0))+126), ones(length(find(mask(chn, :) > 0)), 1) * -.5, 'k', 'LineWidth', 10)
    

    hold off
    box off
    ax = gca;
    axis([-inf 350 -1.1 .7])
    ax.LineWidth = 2.5;
    ax.FontSize = 28;
    ax.FontName = 'Arial';
    xlabel('Time (ms)', 'FontSize', 32, 'FontWeight','normal')
    ylabel('Amplitude (\muV)', 'FontSize', 32, 'FontWeight','normal')
    name_channels =  expected_chanlocs(chn).labels;
    title(name_channels , 'FontSize', 32, 'FontWeight','normal')
        
end


%% data for topoplot - difference across tasks
% find channels where task effect is significant
level2_RM_ANOVA_dir = [pwd, '\level2_RM_ANOVA_HEP_grp_task_effect_50-350ms'];
load([level2_RM_ANOVA_dir, filesep, 'mask_main_effect.mat'])
% average across time window where task effect is significant
time_window = find(sum(mask, 1) > 0);

simpleRT_topoplot = mean(simpleRT_mean(:, 125+time_window), 2);
gng_topoplot = mean(gng_mean(:, 125+time_window), 2);


% channels with significant group difference
sig_chan_number = find(sum(mask, 2) > 0);

title_txt = {'Simple RT' 'GNG'};
       
%         for x = 1:length(sig_chan_number)
%             sig_chans{grp, task}{x} = chanlocs_EEGChanOnly(sig_chan_number(x)).labels;
%         end

figure;
topoplot(ones(59, 1)*.025, expected_chanlocs, 'electrodes', 'off', 'plotchans', sig_chan_number, 'style', 'blank',...
    'plotdisk', 'on',  'hcolor'  , 'none') ; hold on
topoplot(simpleRT_topoplot, expected_chanlocs, 'electrodes', 'off'); 
caxis([-.5 .5]); %c.Axis.FontSize = 16;
cH = colorbar; set(cH,'FontSize',30);
% colorbar('Ticks',[0, 1, 2, 3, 4, 5], 'FontSize', 18, 'FontWeight','normal');
colormap(crameri('vik')); % needs colour maps from http://www.fabiocrameri.ch/colourmaps.php
% title('Young', 'FontSize', 30, 'FontWeight','normal')
% set(get(gca,'title'),'Position',[0,-.65, 0])
%         text(5, 0.4, title_txt{task, grp})


figure;
topoplot(ones(59, 1)*.025, expected_chanlocs, 'electrodes', 'off', 'plotchans', sig_chan_number, 'style', 'blank',...
    'plotdisk', 'on',  'hcolor'  , 'none') ; hold on
topoplot(gng_topoplot, expected_chanlocs, 'electrodes', 'off'); 
caxis([-.5 .5]); %c.Axis.FontSize = 16;
cH = colorbar; set(cH,'FontSize',30);
% colorbar('Ticks',[0, 1, 2, 3, 4, 5], 'FontSize', 18, 'FontWeight','normal');
colormap(crameri('vik')); % needs colour maps from http://www.fabiocrameri.ch/colourmaps.php
% title('Older', 'FontSize', 30, 'FontWeight','normal')
% set(get(gca,'title'),'Position',[0,-.65, 0])


figure;
topoplot(ones(59, 1)*.025, expected_chanlocs, 'electrodes', 'off', 'plotchans', sig_chan_number, 'style', 'blank',...
    'plotdisk', 'on',  'hcolor'  , 'none') ; hold on
topoplot(gng_topoplot-simpleRT_topoplot, expected_chanlocs, 'electrodes', 'off'); 
% caxis([-.4 .4]); %c.Axis.FontSize = 16;
cH = colorbar; set(cH,'FontSize',30);
% colorbar('Ticks',[0, 1, 2, 3, 4, 5], 'FontSize', 18, 'FontWeight','normal');
colormap(crameri('vik')); % needs colour maps from http://www.fabiocrameri.ch/colourmaps.php


%% calculate average amplitude difference across tasks inside clusters where there is a significant task effect
level2_RM_ANOVA_dir = [pwd, '\level2_RM_ANOVA_HEP_grp_task_effect_50-350ms'];
load([level2_RM_ANOVA_dir, filesep, 'mask_main_effect.mat'])

[row, column] = find(mask == 1);

simpleRT_mean = squeeze(mean(HEP_task{1}, 1));
simpleRT_std = squeeze(std(HEP_task{1}, [], 1));

gng_mean = squeeze(mean(HEP_task{2}, 1));
gng_std = squeeze(std(HEP_task{2}, [], 1));

%NOT SURE THIS IS RIGHT!!! CHECK!
avg_ampl_cluster_simpleRT = mean(mean(simpleRT_mean(row, column + 125), 1))
avg_ampl_cluster_gng = mean(mean(gng_mean(row, column + 125), 1))
avg_ampl_cluster_gng-avg_ampl_cluster_simpleRT

HEP_taskdiff = HEP_task{2}-HEP_task{1};
HEP_taskdiff_mean = squeeze(mean(HEP_taskdiff, 1));
HEP_taskdiff_std = squeeze(std(HEP_taskdiff, [], 1));

taskdiff_mean = mean(mean(HEP_taskdiff_mean(row, column + 125), 1))
taskdiff_std = mean(mean(HEP_taskdiff_std(row, column + 125), 1))


% across participants standard deviation
HEP_all_std = squeeze(std(HEP_all, [], 1));
std_ampl_cluster1 = mean(mean(HEP_all_std(row, column + 125), 1))


% effect size use matlab function meanEffectSize


%% average across groups - effect of task
% trimmed means - as used by LIMO toolbox
% TM = limo_trimmed_mean(data,percent)

HEP_task = cell(2, 1);

HEP_task{1} = cat(1, HEP_T_locked{1, 1}, HEP_T_locked{2, 1}); % participants x chan  x time
HEP_task{2} = cat(1, HEP_T_locked{1, 2}, HEP_T_locked{2, 2});

% [t,tmdata,trimci,se,p,tcrit,df]=limo_trimci(data, percent, alpha, nullvalue)
% INPUT: 
%        data (electrodes*frames*trials)
%        percent (percentage of trimming [0-100], default=20)
%        alphav (default=.05)
%        nullvalue (mean of nullhypothesis, default=0)
%
% OUTPUT:
%        t      (2D matrix of t statistics at every electrode and time frame)
%        tmdata (2D matrix of trimmed means)
%        trimci (confidence intervals around the trimmed means)
%        se      standard error
%        p      (p values)
%        tcrit  (1-alpha/2 quantile of the Student's t distribution with
%               adjusted degrees of freedom)
%        df     (degrees of freedom)
%

% rearrange data to have participants in dim 3

for p = 1:size(HEP_task{1}, 1)
    simpleRT_data(:, :, p) =  HEP_task{1}(p, :, :);
end

for p = 1:size(HEP_task{2}, 1)
    gng_data(:, :, p) =  HEP_task{2}(p, :, :);
end


% [t,tmdata,trimci,se,p,tcrit,df]=limo_trimci(data, percent, alpha, nullvalue)

[t,simpleRT_mean,simpleRT_trimci,simpleRT_se,p,tcrit,df]=limo_trimci(simpleRT_data);

[t,gng_mean,gng_trimci,gng_se,p,tcrit,df]=limo_trimci(gng_data);

% simpleRT_trad_mean = squeeze(mean(simpleRT_data, 3));
% 
% 
% figure; plot(simpleRT_trad_mean(34, :)); hold on
% plot(simpleRT_mean(34, :))

% color for each task
clr = [230, 159, 0; 0 114 178; 213, 94, 0]./255;

channels_num = [15 16 23 24 33 34 51];
for chn = channels_num

   
    x_axis = -199:2:400; % in ms
    figure;
    
    % grey background in time window not used for stats
    for idx = -185:5:50
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
    axis([-inf 350 -1.1 1])
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

level2_onesample_dir = [pwd, '\level2_onesample_HEP_grps_together'];
channels_num = [16 34 51];%[7, 25, 43, 57];
for chn = channels_num
    x_axis = -199:2:400; % in ms
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
load([pwd, '\level2_onesample_HEP_grpstogeth_50-350ms\mask.mat'])

[row, column] = find(mask == 1);

avg_ampl_cluster1 = mean(mean(HEP_all_mean(row, column + 125), 1));

% across participants standard deviation
HEP_all_std = squeeze(std(HEP_all, [], 1));
std_ampl_cluster1 = mean(mean(HEP_all_std(row, column + 125), 1))

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
    hold on
    plot([0 0], [-1.2 1], '--k',  'LineWidth', 1) 
    hold on
    
    % horizontal line marking where young group is different from zero
    % find mask for channel chosen   
    level2_onesample_dir = [pwd, '\level2_onesample_HEP_grpstogeth_50-350ms'];
    load([level2_onesample_dir, filesep, 'mask.mat'])
%     plot(x_axis(find((mask(chn, :) == 1))+126), ones(length(find(mask(chn, :) == 1)), 1) * -1, 'm', 'LineWidth', 5)
    plot(x_axis(find(mask(chn, :) > 0) + 126), ones(length(find(mask(chn, :) > 0)), 1) * -.6, 's', 'MarkerSize', 5.5, 'color', 'm', ...
            'MarkerFaceColor', 'm')

% 
%     % horizontal line marking where young group is different from zero
%     % find mask for channel chosen   
%     level2_onesample_dir = [pwd, '\level2_onesample_HEP_separate_grps'];
%     if exist([level2_onesample_dir, filesep, 'young', filesep, 'mask.mat'], 'file')
%         load([level2_onesample_dir, filesep, 'young', filesep, 'mask.mat'])
%         plot(x_axis(find(mask(chn, :) > 0) +126), ones(length(find(mask(chn, :) > 0)), 1) * -.4,'color', [.5 .5 .5], 'LineWidth', 5)
%         hold on
%     end
%     if exist([level2_onesample_dir, filesep, 'older', filesep, 'mask.mat'], 'file')
%         load([level2_onesample_dir, filesep, 'older', filesep, 'mask.mat'])
%         plot(x_axis(find(mask(chn, :) > 0) +126), ones(length(find(mask(chn, :) > 0)), 1) * -.6, 'color', [1 .5 .5], 'LineWidth', 5)
%         hold on
%     end


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


%% plot data using trimmed means
% % [t,tmdata,trimci,se,p,tcrit,df]=limo_trimci(data, percent, alphav, nullvalue)
% % 
% % for s=1:size(HEP{1}, 1)
% %     data_young(:,:, s) = HEP{1}(s,:, :);
% % end
% % 
% % for s=1:size(HEP{2}, 1)
% %     data_older(:,:, s) = HEP{2}(s,:, :);
% % end
% % 
% % [~,young_all_mean,~,young_all_se,~,~,~]=limo_trimci(data_young); % 20% trimmed means
% % 
% % [~,older_all_mean,~,older_all_se,~,~,~]=limo_trimci(data_older); % 20% trimmed means
% % 
% % 
% % channels_num = [16 34 51];%[7, 25, 43, 57];
% % for chn = channels_num
% % 
% %    
% %     x_axis = -199:2:400; % in ms
% %     figure;
% %     
% %     grey background in time window not used for stats
% %     for idx = -185:5:50
% %         plot([idx idx], [-1.15 1], 'color', [.9 .9 .9], 'LineWidth', 5)
% %         hold on
% %     end
% %     
% %     plot(x_axis, young_all_mean(chn, :), 'color', 'k',  'LineWidth',1.5)
% %     hold on
% %     jbfill(x_axis, [young_all_mean(chn, :) + young_all_se(chn, :)], [young_all_mean(chn, :) - young_all_se(chn, :)], 'k', 'k', 1, 0.2)
% %     hold on
% %     plot(x_axis(1:125), young_all_mean(chn, 1:125), 'color', [.6 .6 .6],  'LineWidth',1.5) %plot(light blue)
% %     hold on
% %     plot(x_axis, older_all_mean(chn, :), '--r',  'LineWidth',1.5) 
% %     hold on
% %     jbfill(x_axis, [older_all_mean(chn, :) + older_all_se(chn, :)], [older_all_mean(chn, :) - older_all_se(chn, :)], 'r', 'r', 1, 0.2)
% %     hold on
% %     plot(x_axis(1:125), older_all_mean(chn, 1:125), '--', 'color', [1 .6 .6], 'LineWidth',1.5) 
% %     hold on
% %     plot(x_axis, zeros(length(x_axis), 1), '--k',  'LineWidth', 1) 
% % 
% %     hold off
% %     box off
% %     ax = gca;
% %     axis([-inf 350 -1.2 1])
% %     ax.LineWidth = 2.5;
% %     ax.FontSize = 28;
% %     ax.FontName = 'Arial';
% %     xlabel('Time (ms)', 'FontSize', 32, 'FontWeight','normal')
% %     ylabel('Amplitude (\muV)', 'FontSize', 32, 'FontWeight','normal')
% %     name_channels =  expected_chanlocs(chn).labels;
% %     title(name_channels , 'FontSize', 32, 'FontWeight','normal')
% %         
% % end


%% data for topoplot
% find channels where HEP is significantly different from zero
level2_RM_ANOVA_dir = [pwd filesep 'level2_onesample_HEP_grpstogeth_50-350ms'];
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
caxis([-.5 .5]); %c.Axis.FontSize = 16;
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
caxis([-.5 .5]); %c.Axis.FontSize = 16;
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



