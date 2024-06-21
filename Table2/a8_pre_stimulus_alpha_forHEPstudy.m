% calculation of alpha power in participants included in the
% heartbeat-evoked response - T-locked analyses
% alpha power calculated in the 3 s before cue onset here - 
% G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\pre-stimulus_alpha\fooof_all_channels
% G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\pre-stimulus_alpha\pre_stimulus_power_spectrum.m
% G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\pre-stimulus_alpha\foof_analysis_power_spectrum_thresh1_4pks_all_chans.m
clear
% participants included in original analyses
younger=[4	6	9	10	13	15	16	25	26	28	31	33	34	36	42	44	45	46	50	51	53	54	56	59	62	66	68	72	74	76	78	80	81	82	84	85];
older=[7   8	11	12	14	17	19	20	21	22	23	32	35	37	38	41	43	47	48	49	52	55	57	58	61	63	64	65	67	69	70	71	73	75	77	79	83	86];
participants = {younger, older}; task={'W1', 'D1', 'D2', 'G1', 'G2'};

% participants to exclude that were not included in T-locked HEP analyses
toexclude = [10, 17, 19, 35, 61, 73, 77, 79, 81, 21, 50, 55];

idx2exclude = cell(2, 1);
for grp = 1:2
    for i = toexclude
        idx2exclude{grp} = [idx2exclude{grp}; find(participants{grp} == i)];
    end
end

data_dir = 'G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\pre-stimulus_alpha\fooof_all_channels';
load([data_dir filesep 'oscillations_power_young']); % task x participant x channel x freqband (theta, alpha, beta)
load([data_dir filesep 'oscillations_power_older']);

alpha4HEPstudy{1} = oscillations_power_young(:, setdiff(1:size(oscillations_power_young, 2), idx2exclude{1}), :, 2);
alpha4HEPstudy{2} = oscillations_power_older(:, setdiff(1:size(oscillations_power_older, 2), idx2exclude{2}), :, 2);

%% topoplot average alpha power across participants
% eeglab
load('G:\ProjectAgingNeuromodulation\AuditoryResearch\EEGLAB_analysis\chanlocs.mat');

% for p = 1:size(PowerSpectralDensity_Young, 1)
%     for task = 1:4
%         for channel = 1:size(PowerSpectralDensity_Young, 3)
%             avg_alpha_power_young(p, task, channel) = mean(mean(PowerSpectralDensity_Young{p, task, channel}(5:7, :), 1));
%         end
%     end
% end

% alpha4HEPstudy = task x participant x channel
groups = {'young', 'older'};
for grp = 1:2
    
    % SIMPLE RT
    figure;
    topoplot(squeeze(mean(alpha4HEPstudy{grp}(2, :, :), 2)), chanlocs(1:59))%, 'maplimits', [-8 8]);
    title(['Simple RT - ', groups(grp)]);
    caxis([0 1]); %c.Axis.FontSize = 16;
    colorbar;
    colorbar('Ticks',[0,.5, 1], 'FontSize', 38, 'FontWeight','normal');
    colormap(crameri('imola')); % needs colour maps from http://www.fabiocrameri.ch/colourmaps.php

    % GNG
    figure;
    topoplot(squeeze(mean(alpha4HEPstudy{grp}(3, :, :), 2)), chanlocs(1:59))%, 'maplimits', [-8 8]);
    title(['Go/no-go - ', groups(grp)]);
    caxis([0 1]); %c.Axis.FontSize = 16;
    colorbar;
    colorbar('Ticks',[0,.5, 1], 'FontSize', 38, 'FontWeight','normal');
    colormap(crameri('imola')); % needs colour maps from http://www.fabiocrameri.ch/colourmaps.php
end


%% t-test comparison across tasks

% [pval, t_orig, crit_t, est_alpha, seed_state]=mult_comp_perm_t2(data1,data2)
grp = 1;
[pval, t_orig, crit_t, est_alpha, seed_state]=mult_comp_perm_t2(squeeze(alpha4HEPstudy{grp}(2, :, :)), squeeze(alpha4HEPstudy{grp}(3, :, :)));

find(pval < .05)


[pval, t_orig, crit_t, est_alpha, seed_state]=mult_comp_perm_t2([squeeze(alpha4HEPstudy{1}(2, :, :)); squeeze(alpha4HEPstudy{2}(2, :, :))],...
    [squeeze(alpha4HEPstudy{1}(3, :, :)); squeeze(alpha4HEPstudy{2}(3, :, :))]);

find(pval < .05)


%% Average alpha power in parieto-occipital cluster for repeated measures ANOVA studyign effect of task, effect of group and task x group interaction
% alpha4HEPstudy = task x participant x channel

alpha_simpleRT_young = squeeze(mean(alpha4HEPstudy{1}(2, :, :), 2));
alpha_gng_young = squeeze(mean(alpha4HEPstudy{1}(3, :, :), 2));

alpha_simpleRT_older = squeeze(mean(alpha4HEPstudy{2}(2, :, :), 2));
alpha_gng_older = squeeze(mean(alpha4HEPstudy{2}(3, :, :), 2));

alpha_mean = mean([alpha_gng_older, alpha_gng_young, alpha_simpleRT_older, alpha_simpleRT_young], 2)

[~, idx] = sort(alpha_mean)

chanlocs(idx(end-8:end)).labels

% data for SPSS table

alpha_cluster_simpleRT = squeeze(mean(alpha4HEPstudy{1}(2, :, idx(end-8:end)), 3))';

alpha_cluster_simpleRT = [alpha_cluster_simpleRT; squeeze(mean(alpha4HEPstudy{2}(2, :, idx(end-8:end)), 3))'];

alpha_cluster_gng = squeeze(mean(alpha4HEPstudy{1}(3, :, idx(end-8:end)), 3))';

alpha_cluster_gng = [alpha_cluster_gng; squeeze(mean(alpha4HEPstudy{2}(3, :, idx(end-8:end)), 3))'];

group = [ones(size(alpha4HEPstudy{1}, 2), 1); ones(size(alpha4HEPstudy{2}, 2), 1)*2];

T = array2table([group, alpha_cluster_simpleRT, alpha_cluster_gng],...
    'VariableNames',{'group', 'alpha_simpleRT', 'alpha_gng'});

filename = 'alpha4HEPstudy.xlsx';
writetable(T, filename,'Sheet',1,'Range','A1')

%% plot_2conditions(grp1_cond1, grp1_cond2, grp2_cond1, grp2_cond2, y_label_text, x_label, title_text)

plot_2conditions(alpha_cluster_simpleRT(group == 1), alpha_cluster_gng(group == 1),...
    alpha_cluster_simpleRT(group == 2), alpha_cluster_gng(group == 2), 'Alpha power', {'Simple RT' 'Go/no-go' 'Simple RT' 'Go/no-go'}, '')


function plot_2conditions(grp1_cond1, grp1_cond2, grp2_cond1, grp2_cond2, y_label_text, x_label, title_text)
    figure;
    index_nan1 = find(isnan(grp1_cond1));
    index_nan2 = find(isnan(grp1_cond2));

    grp1_cond1(unique([index_nan1, index_nan2])) = [];
    grp1_cond2(unique([index_nan1, index_nan2])) = [];

    % plot data for condition 1 group 1
    yMean1=nanmean(grp1_cond1);
    y_se = std(grp1_cond1,'omitnan')/sqrt(length(grp1_cond1));
    
        %plot the mean+-SEM box
        %   RECTANGLE('Position',pos) creates a rectangle in 2-D coordinates.
%   Specify pos as a four-element vector of the form [x y w h] in data
%   units. The x and y elements determine the location and the w and h
%   elements determine the size. The function plots into the current axes
%   without clearing existing content from the axes.
    box off
    rectangle('Position',[1-0.3,yMean1-y_se, 0.6, 2*y_se ],'FaceColor',[.7 .7 .7],'EdgeColor', [.7 .7 .7],'LineWidth',0.1);
    hold on
    %plot the mean line    
    plot([1-0.3 1+0.3],[yMean1 yMean1] ,'Color','k','LineWidth',5);
    
    % condition 2 group 1
    yMean2=nanmean(grp1_cond2);
    y_se = std(grp1_cond2,'omitnan')/sqrt(length(grp1_cond2));
    rectangle('Position',[2-0.3,yMean2-y_se, 0.6, 2*y_se ],'FaceColor',[.7 .7 .7],'EdgeColor', [.7 .7 .7],'LineWidth',0.1)
    %plot the mean line
    plot([2-0.3 2+0.3],[yMean2 yMean2] ,'Color','k','LineWidth',5);
    
    for y=1:length(grp1_cond1)
        plot([1+rand*0.2-0.1, 2+rand*0.2-0.1], [grp1_cond1(y), grp1_cond2(y)],'-o', 'color', [.8 .8 .8], ...
            'MarkerEdgeColor','k','MarkerSize',8, 'LineWidth', 1);
        hold on;
    end
    
    % plot data for group 2
    index_nan1 = find(isnan(grp2_cond1));
    index_nan2 = find(isnan(grp2_cond2));

    grp2_cond1(unique([index_nan1, index_nan2])) = [];
    grp2_cond2(unique([index_nan1, index_nan2])) = [];

    % plot data for condition 1
    yMean1=nanmean(grp2_cond1);
    y_se = std(grp2_cond1,'omitnan')/sqrt(length(grp2_cond1));
    
        %plot the mean+-SEM box
        %   RECTANGLE('Position',pos) creates a rectangle in 2-D coordinates.
%   Specify pos as a four-element vector of the form [x y w h] in data
%   units. The x and y elements determine the location and the w and h
%   elements determine the size. The function plots into the current axes
%   without clearing existing content from the axes.
    box off
    rectangle('Position',[4-0.3,yMean1-y_se, 0.6, 2*y_se ],'FaceColor',[.7 .7 .7],'EdgeColor', [.7 .7 .7],'LineWidth',0.1);
    hold on
    %plot the mean line    
    plot([4-0.3 4+0.3],[yMean1 yMean1] ,'Color','k','LineWidth',5);
    
    % condition 2 - group 2
    yMean2=nanmean(grp2_cond2);
    y_se = std(grp2_cond2,'omitnan')/sqrt(length(grp2_cond2));
    rectangle('Position',[5-0.3,yMean2-y_se, 0.6, 2*y_se ],'FaceColor',[.7 .7 .7],'EdgeColor', [.7 .7 .7],'LineWidth',0.1)
    %plot the mean line
    plot([5-0.3 5+0.3],[yMean2 yMean2] ,'Color','k','LineWidth',5);
    
    for y=1:length(grp2_cond1)
        plot([4+rand*0.2-0.1, 5+rand*0.2-0.1], [grp2_cond1(y), grp2_cond2(y)],'-o', 'color', [.8 .8 .8], ...
            'MarkerEdgeColor', 'k','MarkerSize',8, 'LineWidth', 1);
        hold on;
    end

    % axes('XColor','none');
    hold off;
    axis([0 6 -inf inf]);
    xticks([1 2 4 5])
    
    ax = gca;
    ax.LineWidth = 2.5; 
    ax.FontSize = 24;
    ax.FontName = 'Arial';
    ax.Color = 'none';
    ax.XTickLabel= x_label;
    xtickangle(45)
    ax.XAxis.FontSize = 18;

%     ax.XTickLabel= [];

    ylabel(y_label_text, 'FontSize', 26, 'FontWeight','normal')
    title(title_text, 'FontSize', 28, 'FontWeight','normal')

%     x0=10;
%     y0=10;
%     width=400;
%     height=400;
%     set(gcf,'position',[x0,y0,width,height])
end
