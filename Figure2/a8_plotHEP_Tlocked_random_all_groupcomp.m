% plot HEP - all time course
% all participants together and group differences

load HEP_T_locked_random; % calculated in \HEP_grp_task_effect_Tlocked\LIMO_HEPTlocked_grp_task_effect.m
load('D:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\LIMO_stats\expected_chanlocs_both.mat');

% average across tasks
for grp =1:2
    for prt = 1:size(HEP_T_locked_random{grp, 1}, 1)
        clear tmp
        tmp(:, :, 1) = squeeze(HEP_T_locked_random{grp, 1}(prt, :, :));
        tmp(:, :, 2) = squeeze(HEP_T_locked_random{grp, 2}(prt, :, :));
        HEP{grp}(prt, :, :) = mean(tmp, 3); 
    end
end

%% plot HEP both groups together, avg both tasks

HEP_all = cat(1, HEP{1}, HEP{2});
HEP_all_mean = squeeze(mean(HEP_all, 1));
HEP_all_se = squeeze(std(HEP_all, [], 1)/sqrt(size(HEP_all, 1)));

level2_onesample_dir = [pwd, '\level2_onesample_HEP_all_partcts_50-400ms'];
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
    plot(x_axis, zeros(length(x_axis), 1), '--k')
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


%%
young_all_mean = squeeze(mean(HEP{1}, 1));
young_all_se = squeeze(std(HEP{1}, [], 1))/sqrt(size(HEP{1}, 1));

older_all_mean = squeeze(mean(HEP{2}, 1));
older_all_se = squeeze(std(HEP{2}, [], 1))/sqrt(size(HEP{2}, 1));

% % comparison across groups - including time windows with group differences
% level2_RM_ANOVA_dir = [pwd filesep 'level2_RM_ANOVA_HEP_grp_task_effect_50-400ms'];
% load([level2_RM_ANOVA_dir, filesep, 'mask_group_effect.mat'])

channels_num = [16 34 51];%[7, 25, 43, 57];
for chn = channels_num

   
    x_axis = -199:2:400; % in ms
    figure;
    plot(x_axis, young_all_mean(chn, :), 'color', 'k',  'LineWidth',1.5) %plot(light blue)
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
    plot(x_axis, zeros(length(x_axis), 1), '--k')
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
%     plot(x_axis(mask(chn, :) > 0) + 250, ones(length(find(mask(chn, :) > 0)), 1) * -1.15, 's', 'MarkerSize', 4.5, 'color', 'k', ...
%         'MarkerFaceColor', 'k')
    
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
caxis([-1 1]); %c.Axis.FontSize = 16;
cH = colorbar; set(cH,'FontSize',30);
% colorbar('Ticks',[0, 1, 2, 3, 4, 5], 'FontSize', 18, 'FontWeight','normal');
colormap(crameri('vik')); % needs colour maps from http://www.fabiocrameri.ch/colourmaps.php
title('Young', 'FontSize', 30, 'FontWeight','normal')
set(get(gca,'title'),'Position',[0,-.65, 0])
%         text(5, 0.4, title_txt{task, grp})


figure;
topoplot(ones(59, 1)*.025, expected_chanlocs, 'electrodes', 'off', 'plotchans', sig_chan_number, 'style', 'blank',...
    'plotdisk', 'on',  'hcolor'  , 'none') ; hold on
topoplot(older_topoplot, expected_chanlocs, 'electrodes', 'off'); 
caxis([-1 1]); %c.Axis.FontSize = 16;
cH = colorbar; set(cH,'FontSize',30);
% colorbar('Ticks',[0, 1, 2, 3, 4, 5], 'FontSize', 18, 'FontWeight','normal');
colormap(crameri('vik')); % needs colour maps from http://www.fabiocrameri.ch/colourmaps.php
title('Older', 'FontSize', 30, 'FontWeight','normal')
set(get(gca,'title'),'Position',[0,-.65, 0])
     