%% plot betas for each task - coefficients between IBI and HEP amplitude
clear; close all
% create variable with data Y = channel x time x subject

tasks = {'simpleRT', 'gng'};
folder_list = dir(pwd);
load('E:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\LIMO_stats\expected_chanlocs_both.mat');


s = 0; Coefs = [];
for f = 1:length(folder_list)
    if contains(folder_list(f).name, 'AB')
        s = s + 1;
        load([folder_list(f).folder, filesep, folder_list(f).name, filesep, 'Betas.mat'])
        for t = 1:length(tasks)
            Coefs(s, t, :, :) = squeeze(Betas(:, 1:175, t+2));
        end
        subjects{s} = [folder_list(f).name];
    end
end

mean_Coefs = squeeze(mean(Coefs(:, :, :, :), 1));
se_Coefs = squeeze(std(Coefs(:, :, :, :), [], 1))/sqrt(size(Coefs, 1));

%% color for each task
clr = [230, 159, 0; 0 114 178; 213, 94, 0]./255;

% load time windows where the effect of IBI on teh HEP shows an effect of
% task
mask_dir = [pwd, '\level2_RM_ANOVA_IBIvsHEP'];
load([mask_dir, '\mask_main_effect']);

x_axis = 51:2:400;
channels_num = [16, 34, 51];
for chn = channels_num
    figure;
    for t = 1:length(tasks)
        if t == 1
            plot(x_axis, squeeze(mean_Coefs(t, chn, :)), 'color', clr(t+1, :), 'LineWidth',1.5) 
        else
            plot(x_axis, squeeze(mean_Coefs(t, chn, :)), '--', 'color', clr(t+1, :), 'LineWidth',1.5) 
        end
        hold on
        jbfill(x_axis, [squeeze(mean_Coefs(t, chn, :)) + squeeze(se_Coefs(t, chn, :))]',...
            [squeeze(mean_Coefs(t, chn, :)) - squeeze(se_Coefs(t, chn, :))]',...
            clr(t+1, :), clr(t+1, :), 1, .2)
        hold on
    end
    plot(x_axis, zeros(length(x_axis), 1), '--k', 'LineWidth',1) 
    
    mask_idx = find(mask(chn, :) == 1);
    plot(x_axis(mask_idx), ones(length(mask_idx), 1)  * -.2, 's', 'MarkerSize', 4.5, 'color', 'k', ...
            'MarkerFaceColor', 'k')
    
    axis([51 400 -.3 .6])
    box off
    ax = gca;
    ax.LineWidth = 2.5;
    ax.FontSize = 28;
    ax.FontName = 'Arial';
    xlabel('Time (ms)', 'FontSize', 32, 'FontWeight','normal')
    ylabel('Coefficients', 'FontSize', 32, 'FontWeight','normal')
    name_channels =  expected_chanlocs(chn).labels;
    title(name_channels, 'FontSize', 32, 'FontWeight','normal')
%     legend('simpleRT', 'gng')
end

%% plot coefficients young vs older
older = [11 12 14 20 21 22 23 32 37 38 41 43 47 48 49 52 55 57 58 63 64 65 67 69 7 70 71 75 8 83 86];
young = [13 15 16 25 26 28 31 33 34 36 4 42 44 45 46 50 51 53 54 56 59 6 62 66 68 72 74 76 78 80 82 84 85 9];
Coefs_young = []; Coefs_older = [];
for t = 1:length(tasks)
    yng = 0; old = 0; 
    for f = 1:length(folder_list)
        if contains(folder_list(f).name, 'AB')
            load([folder_list(f).folder, filesep, folder_list(f).name, filesep, 'Betas.mat'])
            if ismember(str2num(folder_list(f).name(3:end)), young)
                yng = yng + 1;
                Coefs_young( yng, t, :, :) = squeeze(Betas(:, 1:175, t+2));
            else
                old = old +1;
                Coefs_older( old, t, :, :) = squeeze(Betas(:, 1:175, t+2));
            end
        end
    end
end


mean_Coefs_young = squeeze(mean(Coefs_young(:, :, :, :), 1));
se_Coefs_young = squeeze(std(Coefs_young(:, :, :, :), [], 1))/sqrt(size(Coefs_young, 1));

mean_Coefs_older = squeeze(mean(Coefs_older(:, :, :, :), 1));
se_Coefs_older = squeeze(std(Coefs_older(:, :, :, :), [], 1))/sqrt(size(Coefs_older, 1));

%% color for each task
clr = [230, 159, 0; 0 114 178; 213, 94, 0]./255;

% load time windows where the effect of IBI on teh HEP shows an effect of
% task
% mask_dir = [pwd, '\level2_RM_ANOVA_IBIvsHEP'];
% load([mask_dir, '\mask_main_effect']);

x_axis = 51:2:400;
channels_num = [16, 34, 51];
for t = 1:length(tasks)
    for chn = channels_num
        figure;
        plot(x_axis, squeeze(mean_Coefs_young(t, chn, :)), 'color', 'k', 'LineWidth',1.5) 
        hold on
        jbfill(x_axis, [squeeze(mean_Coefs_young(t, chn, :)) + squeeze(se_Coefs_young(t, chn, :))]',...
            [squeeze(mean_Coefs_young(t, chn, :)) - squeeze(se_Coefs_young(t, chn, :))]',...
            'k', 'k', 1, .2)
        hold on
        plot(x_axis, squeeze(mean_Coefs_older(t, chn, :)), '--', 'color', 'r', 'LineWidth',1.5) 
        hold on
        jbfill(x_axis, [squeeze(mean_Coefs_older(t, chn, :)) + squeeze(se_Coefs_older(t, chn, :))]',...
            [squeeze(mean_Coefs_older(t, chn, :)) - squeeze(se_Coefs_older(t, chn, :))]',...
            'r', 'r', 1, .2)
        hold on
        plot(x_axis, zeros(length(x_axis), 1), '--k', 'LineWidth',1) 

%         mask_idx = find(mask(chn, :) == 1);
%         plot(x_axis(mask_idx), ones(length(mask_idx), 1)  * -.2, 's', 'MarkerSize', 4.5, 'color', 'k', ...
%                 'MarkerFaceColor', 'k')

        axis([51 400 -inf inf])
        box off
        ax = gca;
        ax.LineWidth = 2.5;
        ax.FontSize = 28;
        ax.FontName = 'Arial';
        xlabel('Time (ms)', 'FontSize', 32, 'FontWeight','normal')
        ylabel('Coefficients', 'FontSize', 32, 'FontWeight','normal')
        name_channels =  expected_chanlocs(chn).labels;
        title(name_channels, 'FontSize', 32, 'FontWeight','normal')
    %     legend('simpleRT', 'gng')
    end
end


%% young both tasks

% load time windows where the effect of IBI on teh HEP shows an effect of
% task
% mask_dir = [pwd, '\level2_RM_ANOVA_IBIvsHEP'];
% load([mask_dir, '\mask_main_effect']);

x_axis = 51:2:400;
channels_num = [16, 34, 51];

for chn = channels_num
    figure;
    t=1;
    plot(x_axis, squeeze(mean_Coefs_young(t, chn, :)), 'color', 'k', 'LineWidth',1.5) 
    hold on
    jbfill(x_axis, [squeeze(mean_Coefs_young(t, chn, :)) + squeeze(se_Coefs_young(t, chn, :))]',...
        [squeeze(mean_Coefs_young(t, chn, :)) - squeeze(se_Coefs_young(t, chn, :))]',...
        'k', 'k', 1, .2)
    hold on
    t=2;
    plot(x_axis, squeeze(mean_Coefs_young(t, chn, :)), '--','color', 'b', 'LineWidth',1.5) 
    hold on
    jbfill(x_axis, [squeeze(mean_Coefs_young(t, chn, :)) + squeeze(se_Coefs_young(t, chn, :))]',...
        [squeeze(mean_Coefs_young(t, chn, :)) - squeeze(se_Coefs_young(t, chn, :))]',...
        'b', 'b', 1, .2)
    hold on
    plot(x_axis, zeros(length(x_axis), 1), '--k', 'LineWidth',1) 

%         mask_idx = find(mask(chn, :) == 1);
%         plot(x_axis(mask_idx), ones(length(mask_idx), 1)  * -.2, 's', 'MarkerSize', 4.5, 'color', 'k', ...
%                 'MarkerFaceColor', 'k')

    axis([51 400 -inf inf])
    box off
    ax = gca;
    ax.LineWidth = 2.5;
    ax.FontSize = 28;
    ax.FontName = 'Arial';
    xlabel('Time (ms)', 'FontSize', 32, 'FontWeight','normal')
    ylabel('Coefficients', 'FontSize', 32, 'FontWeight','normal')
    name_channels =  expected_chanlocs(chn).labels;
    title(name_channels, 'FontSize', 32, 'FontWeight','normal')
%     legend('simpleRT', 'gng')
end

%% older both tasks
for chn = channels_num
    figure;
    t=1;
    plot(x_axis, squeeze(mean_Coefs_older(t, chn, :)), 'color', 'k', 'LineWidth',1.5) 
    hold on
    jbfill(x_axis, [squeeze(mean_Coefs_older(t, chn, :)) + squeeze(se_Coefs_older(t, chn, :))]',...
        [squeeze(mean_Coefs_older(t, chn, :)) - squeeze(se_Coefs_older(t, chn, :))]',...
        'k', 'k', 1, .2)
    hold on
    t=2;
    plot(x_axis, squeeze(mean_Coefs_older(t, chn, :)), '--','color', 'b', 'LineWidth',1.5) 
    hold on
    jbfill(x_axis, [squeeze(mean_Coefs_older(t, chn, :)) + squeeze(se_Coefs_older(t, chn, :))]',...
        [squeeze(mean_Coefs_older(t, chn, :)) - squeeze(se_Coefs_older(t, chn, :))]',...
        'b', 'b', 1, .2)
    hold on
    plot(x_axis, zeros(length(x_axis), 1), '--k', 'LineWidth',1) 

%         mask_idx = find(mask(chn, :) == 1);
%         plot(x_axis(mask_idx), ones(length(mask_idx), 1)  * -.2, 's', 'MarkerSize', 4.5, 'color', 'k', ...
%                 'MarkerFaceColor', 'k')

    axis([51 400 -inf inf])
    box off
    ax = gca;
    ax.LineWidth = 2.5;
    ax.FontSize = 28;
    ax.FontName = 'Arial';
    xlabel('Time (ms)', 'FontSize', 32, 'FontWeight','normal')
    ylabel('Coefficients', 'FontSize', 32, 'FontWeight','normal')
    name_channels =  expected_chanlocs(chn).labels;
    title(name_channels, 'FontSize', 32, 'FontWeight','normal')
%     legend('simpleRT', 'gng')
end

%% plot coefficients for participants with high/low pupil dilation responses

pupil_dir = 'E:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\LIMO_stats\HEP_TLockedvsIBI\pupil_dilation_amplitude\';

load([pupil_dir, 'BlnDivided_Part_PeakAmpLat_Median_D_CueLocked_YoungGrp']);
load([pupil_dir, 'BlnDivided_Part_PeakAmpLat_Median_G_CueLocked_YoungGrp']);
load([pupil_dir, 'BlnDivided_Part_PeakAmpLat_Median_D_CueLocked_OlderGrp']);
load([pupil_dir, 'BlnDivided_Part_PeakAmpLat_Median_G_CueLocked_OlderGrp']);


% sort participants according to pupil dilation amplitude
simpleRT_pupil = [BlnDivided_Part_PeakAmpLat_Median_D_CueLocked_YoungGrp;...
    BlnDivided_Part_PeakAmpLat_Median_D_CueLocked_OlderGrp];

gng_pupil = [BlnDivided_Part_PeakAmpLat_Median_G_CueLocked_YoungGrp;...
    BlnDivided_Part_PeakAmpLat_Median_G_CueLocked_OlderGrp];

for s = 1:length(subjects)
    subj_numb(s) = str2num(subjects{s}(3:end));
end

% delete subjects in pupil but not included in HEP analysis
subj_delete = [];
for s = 1:size(simpleRT_pupil, 1)
    if ~ismember(simpleRT_pupil(s, 1), subj_numb)
        subj_delete = [subj_delete; s simpleRT_pupil(s, 1)];
    end
end

simpleRT_pupil(subj_delete(:, 1), :) = [];
gng_pupil(subj_delete(:, 1), :) = [];

simpleRT_pupil_sorted = sortrows(simpleRT_pupil, 2);
gng_pupil_sorted = sortrows(gng_pupil, 2);

coefs_lowpupil_simpleRT = NaN(floor(size(simpleRT_pupil, 1)/2), size(Coefs, 3), size(Coefs, 4));

for s = 1:floor(size(simpleRT_pupil_sorted, 1)/2)
    coefs_lowpupil_simpleRT(s, :, :) = Coefs(subj_numb == simpleRT_pupil_sorted(s, 1), 1, :, :);
end

coefs_highpupil_simpleRT = NaN(floor(size(simpleRT_pupil, 1)/2), size(Coefs, 3), size(Coefs, 4));
idx = 0;
for s = floor(size(simpleRT_pupil_sorted, 1)/2)+1:size(simpleRT_pupil_sorted, 1)
    idx = idx + 1;
    coefs_highpupil_simpleRT(idx, :, :) = Coefs(subj_numb == simpleRT_pupil_sorted(s, 1), 1, :, :);
end


coefs_lowpupil_gng = NaN(floor(size(gng_pupil, 1)/2), size(Coefs, 3), size(Coefs, 4));

for s = 1:floor(size(gng_pupil_sorted, 1)/2)
    coefs_lowpupil_gng(s, :, :) = Coefs(subj_numb == gng_pupil_sorted(s, 1), 1, :, :);
end

coefs_highpupil_gng = NaN(floor(size(gng_pupil, 1)/2), size(Coefs, 3), size(Coefs, 4));
idx = 0;
for s = floor(size(gng_pupil_sorted, 1)/2)+1:size(gng_pupil_sorted, 1)
    idx = idx + 1;
    coefs_highpupil_gng(idx, :, :) = Coefs(subj_numb == gng_pupil_sorted(s, 1), 1, :, :);
end

% color for each task
clr = [230, 159, 0; 0 114 178; 213, 94, 0]./255;

% load time windows where the effect of IBI on teh HEP shows an effect of
% task
% mask_dir = [pwd, '\level2_RM_ANOVA_IBIvsHEP'];
% load([mask_dir, '\mask_main_effect']);

x_axis = 51:2:400;
channels_num = [7, 16, 34, 51, 57];
for chn = channels_num
    figure;
        plot(x_axis, squeeze(mean(coefs_highpupil_simpleRT(:, chn, :), 1)), 'color', clr(2, :), 'LineWidth',1.5) 
        hold on
        plot(x_axis, squeeze(mean(coefs_lowpupil_simpleRT(:, chn, :), 1)), '--', 'color', clr(3, :), 'LineWidth',1.5) 

%         hold on
%         jbfill(x_axis, [squeeze(mean_Coefs(t, chn, :)) + squeeze(se_Coefs(t, chn, :))]',...
%             [squeeze(mean_Coefs(t, chn, :)) - squeeze(se_Coefs(t, chn, :))]',...
%             clr(t+1, :), clr(t+1, :), 1, .2)
        hold on
% end
    plot(x_axis, zeros(length(x_axis), 1), '--k', 'LineWidth',1) 
    
%     mask_idx = find(mask(chn, :) == 1);
%     plot(x_axis(mask_idx), ones(length(mask_idx), 1)  * -.2, 's', 'MarkerSize', 4.5, 'color', 'k', ...
%             'MarkerFaceColor', 'k')
    
    axis([51 400 -.3 .6])
    box off
    ax = gca;
    ax.LineWidth = 2.5;
    ax.FontSize = 28;
    ax.FontName = 'Arial';
    xlabel('Time (ms)', 'FontSize', 32, 'FontWeight','normal')
    ylabel('Coefficients', 'FontSize', 32, 'FontWeight','normal')
    name_channels =  expected_chanlocs(chn).labels;
    title(name_channels, 'FontSize', 32, 'FontWeight','normal')
%     legend('simpleRT', 'gng')
end
%%
x_axis = 51:2:400;
channels_num = [7, 16, 34, 51, 57];
for chn = channels_num
    figure;
        plot(x_axis, squeeze(mean(coefs_highpupil_gng(:, chn, :), 1)), 'color', clr(2, :), 'LineWidth',1.5) 
        hold on
        plot(x_axis, squeeze(mean(coefs_lowpupil_gng(:, chn, :), 1)), '--', 'color', clr(3, :), 'LineWidth',1.5) 

%         hold on
%         jbfill(x_axis, [squeeze(mean_Coefs(t, chn, :)) + squeeze(se_Coefs(t, chn, :))]',...
%             [squeeze(mean_Coefs(t, chn, :)) - squeeze(se_Coefs(t, chn, :))]',...
%             clr(t+1, :), clr(t+1, :), 1, .2)
        hold on
% end
    plot(x_axis, zeros(length(x_axis), 1), '--k', 'LineWidth',1) 
    
%     mask_idx = find(mask(chn, :) == 1);
%     plot(x_axis(mask_idx), ones(length(mask_idx), 1)  * -.2, 's', 'MarkerSize', 4.5, 'color', 'k', ...
%             'MarkerFaceColor', 'k')
    
    axis([51 400 -.3 .6])
    box off
    ax = gca;
    ax.LineWidth = 2.5;
    ax.FontSize = 28;
    ax.FontName = 'Arial';
    xlabel('Time (ms)', 'FontSize', 32, 'FontWeight','normal')
    ylabel('Coefficients', 'FontSize', 32, 'FontWeight','normal')
    name_channels =  expected_chanlocs(chn).labels;
    title(name_channels, 'FontSize', 32, 'FontWeight','normal')
%     legend('gng', 'gng')
end





