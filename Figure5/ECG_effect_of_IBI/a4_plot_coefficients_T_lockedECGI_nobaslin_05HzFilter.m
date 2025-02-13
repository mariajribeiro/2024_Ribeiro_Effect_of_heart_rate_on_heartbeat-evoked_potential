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
            Coefs(s, t, :) = squeeze(Betas(:, :, t+2));
        end
        subjects{s} = [folder_list(f).name];
    end
end

mean_Coefs = squeeze(mean(Coefs(:, :, 1:150), 1));
se_Coefs = squeeze(std(Coefs(:, :, 1:150), [], 1))/sqrt(size(Coefs, 1));

%% color for each task
clr = [230, 159, 0; 0 114 178; 213, 94, 0]./255;

% load time windows where the effect of IBI on teh HEP is significantly
% different from zero for each task condition
% simple RT
mask_dir = [pwd, '\level2_onesample_analysis_50-350ms\simpleRT'];
load([mask_dir, '\mask']);
mask_tasks{1} = mask;
% gng
mask_dir = [pwd, '\level2_onesample_analysis_50-350ms\gng'];
load([mask_dir, '\mask']);
mask_tasks{2} = mask;

x_axis = 51:2:350;
channels_num = 1;
for chn = channels_num
    figure;
    for t = 1:length(tasks)
        if t == 1
            plot(x_axis, squeeze(mean_Coefs(t, :)), 'color', clr(t+1, :), 'LineWidth',1.5) 
        else
            plot(x_axis, squeeze(mean_Coefs(t, :)), '--', 'color', clr(t+1, :), 'LineWidth',1.5) 
        end
        hold on
        jbfill(x_axis, [squeeze(mean_Coefs(t, :)) + squeeze(se_Coefs(t, :))],...
            [squeeze(mean_Coefs(t, :)) - squeeze(se_Coefs(t, :))],...
            clr(t+1, :), clr(t+1, :), 1, .2)
        hold on
        mask_idx = find(mask_tasks{t} == 1);
        if t == 1
            location = -1.5;
        elseif t == 2
            location = -2;
        end
        plot(x_axis(mask_idx), ones(length(mask_idx), 1)  * location, 's', 'MarkerSize', 4.5, 'color', clr(t+1, :), ...
            'MarkerFaceColor', clr(t+1, :))
        hold on
    end
    plot(x_axis, zeros(length(x_axis), 1), '--k', 'LineWidth', 1.5)
    
    axis([51 350 -inf inf])
    box off
    ax = gca;
    ax.LineWidth = 2.5;
    ax.FontSize = 28;
    ax.FontName = 'Arial';
    xlabel('Time (ms)', 'FontSize', 32, 'FontWeight','normal')
    ylabel('Coefficients', 'FontSize', 32, 'FontWeight','normal')
    title('ECG', 'FontSize', 32, 'FontWeight','normal')
%     legend('simpleRT', 'gng')
end


%% plot coefficients with trimmed means from one sample test output

% color for each task
clr = [230, 159, 0; 0 114 178; 213, 94, 0]./255;

tasks = {'simpleRT', 'gng'};

% load time windows where the effect of IBI on teh HEP is significantly
% different from zero for each task condition
% simple RT
mask_dir = [pwd, '\level2_onesample_analysis_50-350ms\simpleRT'];
load([mask_dir, '\mask']);
load([mask_dir, '\one_sample_ttest_parameter_1']);
one_sample_tasks{1} = one_sample;
% gng
mask_dir = [pwd, '\level2_onesample_analysis_50-350ms\gng'];
load([mask_dir, '\mask']);
mask_tasks{2} = mask;
load([mask_dir, '\one_sample_ttest_parameter_1']);
one_sample_tasks{2} = one_sample;

x_axis = 51:2:350;
channels_num = 1;
for chn = channels_num
    figure;
    for t = 1:length(tasks)
        if t == 1
            plot(x_axis, squeeze(one_sample_tasks{t}(1, :, 1)), 'color', clr(t+1, :), 'LineWidth',1.5) 
        else
            plot(x_axis, squeeze(one_sample_tasks{t}(1, :, 1)), '--','color', clr(t+1, :), 'LineWidth',1.5) 
        end
        hold on
        jbfill(x_axis, squeeze(one_sample_tasks{t}(1, :, 1)) + squeeze(one_sample_tasks{t}(1, :, 2)),...
            squeeze(one_sample_tasks{t}(1, :, 1)) - squeeze(one_sample_tasks{t}(1, :, 2)),...
            clr(t+1, :), clr(t+1, :), 1, .2)
        hold on
        mask_idx = find(mask_tasks{t} == 1);
        if t == 1
            location = -1;
        elseif t == 2
            location = -1.5;
        end
        plot(x_axis(mask_idx), ones(length(mask_idx), 1)  * location, 's', 'MarkerSize', 4.5, 'color', clr(t+1, :), ...
            'MarkerFaceColor', clr(t+1, :))
        hold on
    end
    plot(x_axis, zeros(length(x_axis), 1), '--k', 'LineWidth',1) 
    
    axis([51 350 -inf inf])
    ylim padded
    box off
    ax = gca;
    ax.LineWidth = 2.5;
    ax.FontSize = 28;
    ax.FontName = 'Arial';
    xlabel('Time (ms)', 'FontSize', 32, 'FontWeight','normal')
    ylabel('Coefficients', 'FontSize', 32, 'FontWeight','normal')
    title('ECG', 'FontSize', 32, 'FontWeight','normal')
%     legend('simpleRT', 'gng')
end


%% plot coefficients young vs older - average of both tasks
folder_list = dir(pwd);
load('E:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\LIMO_stats\expected_chanlocs_both.mat');
tasks = {'simpleRT', 'gng'};
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
                Coefs_young( yng, t, :, :) = squeeze(Betas(:, :, t+2));
            else
                old = old +1;
                Coefs_older( old, t, :, :) = squeeze(Betas(:, :, t+2));
            end
        end
    end
end


mean_Coefs_young = squeeze(mean(mean(Coefs_young(:, :, :, :), 2), 1));
se_Coefs_young = squeeze(std(mean(Coefs_young(:, :, :, :), 2), [], 1))/sqrt(size(Coefs_young, 1));

mean_Coefs_older = squeeze(mean(mean(Coefs_older(:, :, :, :), 2), 1));
se_Coefs_older = squeeze(std(mean(Coefs_older(:, :, :, :), 2), [], 1))/sqrt(size(Coefs_older, 1));

% color for each task
clr = [230, 159, 0; 0 114 178; 213, 94, 0]./255;


x_axis = 51:2:400;

figure;
plot(x_axis, mean_Coefs_young, 'color', 'k', 'LineWidth',1.5) 
hold on
jbfill(x_axis, [squeeze(mean_Coefs_young) + squeeze(se_Coefs_young)]',...
    [squeeze(mean_Coefs_young) - squeeze(se_Coefs_young)]',...
    'k', 'k', 1, .2)
hold on
plot(x_axis, squeeze(mean_Coefs_older), '--', 'color', 'r', 'LineWidth',1.5) 
hold on
jbfill(x_axis, [squeeze(mean_Coefs_older) + squeeze(se_Coefs_older)]',...
    [squeeze(mean_Coefs_older) - squeeze(se_Coefs_older)]',...
    'r', 'r', 1, .2)
hold on
plot(x_axis, zeros(length(x_axis), 1), '--k', 'LineWidth',1) 


    axis([51 350 -inf inf])%axis([51 400 -.3 .2])
box off
ax = gca;
ax.LineWidth = 2.5;
ax.FontSize = 28;
ax.FontName = 'Arial';
xlabel('Time (ms)', 'FontSize', 32, 'FontWeight','normal')
ylabel('Coefficients', 'FontSize', 32, 'FontWeight','normal')
% name_channels =  expected_chanlocs(chn).labels;
% title(name_channels, 'FontSize', 32, 'FontWeight','normal')
%     legend('simpleRT', 'gng')

 

