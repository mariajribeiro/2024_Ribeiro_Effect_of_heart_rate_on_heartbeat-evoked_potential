%% plot coefficients for participants with high/low pupil dilation responses
% plot betas for each task - coefficients between IBI and HEP amplitude
clear; close all
% create variable with data Y = channel x time x subject

% amplitude of pupil phasic responses
pupil_dir = 'D:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\LIMO_stats\HEP_TLockedvsIBI\pupil_dilation_amplitude\';

load([pupil_dir, 'BlnDivided_Part_PeakAmpLat_Median_D_CueLocked_YoungGrp']);
load([pupil_dir, 'BlnDivided_Part_PeakAmpLat_Median_G_CueLocked_YoungGrp']);
load([pupil_dir, 'BlnDivided_Part_PeakAmpLat_Median_D_CueLocked_OlderGrp']);
load([pupil_dir, 'BlnDivided_Part_PeakAmpLat_Median_G_CueLocked_OlderGrp']);


pupil_amplitude{1} = [BlnDivided_Part_PeakAmpLat_Median_D_CueLocked_YoungGrp(:, 1:2); ...
    BlnDivided_Part_PeakAmpLat_Median_D_CueLocked_OlderGrp(:, 1:2)];

pupil_amplitude{2} = [BlnDivided_Part_PeakAmpLat_Median_G_CueLocked_YoungGrp(:, 1:2); ...
    BlnDivided_Part_PeakAmpLat_Median_G_CueLocked_OlderGrp(:, 1:2)];

pupil_amplitude_diff = pupil_amplitude{2}(:, 2)-pupil_amplitude{1}(:, 2);


tasks = {'simpleRT', 'gng'};
folder_list = dir(pwd);

s = 0; CoefsDiff = []; pupil = [];
for f = 1:length(folder_list)
    if contains(folder_list(f).name, 'AB')
        % check if participant has pupil data
        subj_idx = find(pupil_amplitude{1}(:, 1) == str2num(folder_list(f).name(3:end)));
        if ~isempty(subj_idx)
            s = s + 1;
            pupil(s, 1) = pupil_amplitude_diff(subj_idx, 1);
            load([folder_list(f).folder, filesep, folder_list(f).name, filesep, 'Betas.mat'])
            % diff in regression coefficients across tasks - gng-simpleRT
            CoefsDiff(s, :, :) = Betas(:, 1:175, 4)-Betas(:, 1:175, 3); %channels, time, subjects)
            subjects{s} = [folder_list(f).name];
        end
    end
end


[pupil_sorted, srt_idx] = sortrows(pupil);

coefs_lowdiffpupil = CoefsDiff(srt_idx(1:floor(length(pupil)/2)), :, :);
coefs_highdiffpupil = CoefsDiff(srt_idx(floor(length(pupil)/2)+1:end), :, :);


%% color for each task
load('D:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\LIMO_stats\expected_chanlocs_both.mat');

clr = [230, 159, 0; 0 114 178; 213, 94, 0]./255;

% load time windows where the effect of IBI on the HEP shows an effect of
% task
mask_dir = [pwd, '\level2_regression_pupil_diff'];
load([mask_dir, '\mask_covariate']);



x_axis = 51:2:400;
channels_num = [9 18 27 36 45]% [7, 16, 34, 51, 57];
for chn = channels_num
    
    mean_high = squeeze(mean(coefs_highdiffpupil(:, chn, :), 1));
    se_high = squeeze(std(coefs_highdiffpupil(:, chn, :), [], 1))/sqrt(size(coefs_highdiffpupil, 1));

    mean_low = squeeze(mean(coefs_lowdiffpupil(:, chn, :), 1));
    se_low = squeeze(std(coefs_lowdiffpupil(:, chn, :), [], 1))/sqrt(size(coefs_lowdiffpupil, 1));

       
    figure;
        plot(x_axis, mean_high', 'color', [0.4940 0.1840 0.5560], 'LineWidth',1.5) 
        hold on
        jbfill(x_axis, [mean_high + se_high]',[mean_high - se_high]', [0.4940 0.1840 0.5560], [0.4940 0.1840 0.5560], 1, .2)
        hold on
        plot(x_axis, mean_low', '--', 'color', [0.4660 0.6740 0.1880], 'LineWidth',1.5) 
        hold on
        jbfill(x_axis, [mean_low + se_low]',[mean_low - se_low]', [0.4660 0.6740 0.1880], [0.4660 0.6740 0.1880], 1, .2)
        hold on
        plot(x_axis, zeros(length(x_axis), 1), '--k', 'LineWidth',1) 
    
        mask_idx = find(mask(chn, :) > 0);
        plot(x_axis(mask_idx), ones(length(mask_idx), 1)  * -.25, 's', 'MarkerSize', 4.5, 'color', 'k', ...
                'MarkerFaceColor', 'k')
    
    axis([51 400 -.3 inf])
    box off
    ax = gca;
    ax.LineWidth = 2.5;
    ax.FontSize = 28;
    ax.FontName = 'Arial';
    xlabel('Time (ms)', 'FontSize', 32, 'FontWeight','normal')
    ylabel('\Delta coefficients', 'FontSize', 32, 'FontWeight','normal')
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