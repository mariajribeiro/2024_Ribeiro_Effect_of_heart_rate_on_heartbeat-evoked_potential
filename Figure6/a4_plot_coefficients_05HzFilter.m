%% plot betas for each task - coefficients between IBI and HEP amplitude
clear; close all
% create variable with data Y = channel x time x subject

tasks = {'simpleRT', 'gng'};
folder_list = dir(pwd);
load('D:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\LIMO_stats\expected_chanlocs_both.mat');


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

% % load time windows where the effect of IBI on teh HEP shows an effect of
% % task
% mask_dir = [pwd, '\level2_RM_ANOVA_IBIvsHEP'];
% load([mask_dir, '\mask_main_effect']);

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
    
%     mask_idx = find(mask(chn, :) == 1);
%     plot(x_axis(mask_idx), ones(length(mask_idx), 1)  * -.2, 's', 'MarkerSize', 4.5, 'color', 'k', ...
%             'MarkerFaceColor', 'k')
    
    axis([51 400 -.3 .2])
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
