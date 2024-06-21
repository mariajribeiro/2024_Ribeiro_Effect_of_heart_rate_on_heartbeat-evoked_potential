
% plot HEP inchannels FCz and CPz for participants wuth high and low heart
% rate
clear; close all
tasks = {'simpleRT', 'gng'}; group = {'young', 'older'};
load('D:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\LIMO_stats\expected_chanlocs_both.mat');

% determine participants with high and low heart rate
% load HR data calculated in D:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\hep_effect_of_heart_rate.m
load('D:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\HR_avg_std_min_max_range.mat');

folder_list = dir(pwd);
for grp = 1:2
    for t = 1:length(tasks)

        s = 0; Y = []; HR = []; Yavg = []; HRavg = [];
        for f = 1:length(folder_list)
            if contains(folder_list(f).name, 'AB')
                % check if participant has HR data
                subj_idx = find(HR_avg_std_min_max_range{grp, 1} == str2num(folder_list(f).name(3:end)));
                if ~isempty(subj_idx)
                    s = s + 1;
                    HR(s, 1) = HR_avg_std_min_max_range{grp, 2}(subj_idx, t, 1);
                    HRavg(s, 1) = mean([HR_avg_std_min_max_range{grp, 2}(subj_idx, 2, 1), ...
                        HR_avg_std_min_max_range{grp, 2}(subj_idx, 3, 1)]);
                    load([folder_list(f).folder, filesep, folder_list(f).name, filesep, 'HEP_T_locked_data.mat'])
                    load([folder_list(f).folder, filesep, folder_list(f).name, filesep, 'categorical_variable.mat'])
                    HEP_task = HEP_T_locked_data(:, :, categorical_variable == t);
                    Y(:, :, s) = mean(HEP_task, 3); %channels, time, subjects)
                    Yavg(:, :, s) = mean(HEP_T_locked_data, 3);
                    subjects{s} = [folder_list(f).name];
                end
            end
        end
        HEP_highHR{grp, t} = Y(:, :, HR > median(HR));
        HEP_lowHR{grp, t} = Y(:, :, HR <= median(HR));
        HEPavg_highHR{grp} = Yavg(:, :, HRavg > median(HRavg));
        HEPavg_lowHR{grp} = Yavg(:, :, HRavg <= median(HRavg));
    end

    
end


%% average across tasks
groups = {'young', 'older'};
dir_HR = 'D:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\LIMO_stats\HEP_TLockedvsIBI_nobaseline_05HzFilter\level2_regression_HR_avgtasks';
for grp = 1:2
    mean_highHR = mean(HEPavg_highHR{grp}, 3);
    se_highHR = std(HEPavg_highHR{grp}, [], 3)/sqrt(size(HEPavg_highHR{grp}, 3));

    mean_lowHR = mean(HEPavg_lowHR{grp}, 3);
    se_lowHR = std(HEPavg_lowHR{grp}, [], 3)/sqrt(size(HEPavg_lowHR{grp}, 3));

    channels_num = [7 16 34 51];%[7, 25, 43, 57];
    for chn = channels_num


        x_axis = -199:2:400; % in ms
        figure;
        plot(x_axis,mean_highHR(chn, :), 'color', [0 158 115]./255,  'LineWidth',1.5)
        hold on
        jbfill(x_axis, [mean_highHR(chn, :) + se_highHR(chn, :)], [mean_highHR(chn, :) - se_highHR(chn, :)], [0 158 115]./255, [0 158 115]./255, 1, 0.2)
        hold on
        plot(x_axis, mean_lowHR(chn, :), '--', 'color', [204 121 167]./255,  'LineWidth',1.5) 
        hold on
        jbfill(x_axis, [mean_lowHR(chn, :) + se_lowHR(chn, :)], [mean_lowHR(chn, :) - se_lowHR(chn, :)],[204 121 167]./255, [204 121 167]./255, 1, 0.2)
        hold on
        plot(x_axis, zeros(length(x_axis), 1), '--k',  'LineWidth', 1) 
        hold on

        % horizontal line marking the cluster showing an effect of heart
        % rate
        mask_dir = [dir_HR, filesep, groups{grp}];
        load([mask_dir, '\mask_covariate.mat'])
        plot(x_axis(find(mask(chn, :) > 0) + 125), ones(length(find(mask(chn, :) == 1)), 1) * -1.25, 'k', 'LineWidth', 4.5)
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
        title([groups{grp}, ' ', name_channels] , 'FontSize', 32, 'FontWeight','normal')

    end
end


%% gng task
groups = {'young', 'older'};
dir_HR = 'D:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\LIMO_stats\HEP_TLockedvsIBI_nobaseline_05HzFilter\level2_regression_HR';
for grp = 1:2
    mean_highHR = mean(HEP_highHR{grp, 2}, 3);
    se_highHR = std(HEP_highHR{grp, 2}, [], 3)/sqrt(size(HEP_highHR{grp, 2}, 3));

    mean_lowHR = mean(HEP_lowHR{grp, 2}, 3);
    se_lowHR = std(HEP_lowHR{grp, 2}, [], 3)/sqrt(size(HEP_lowHR{grp, 2}, 3));

    channels_num = [7 16 34 51];%[7, 25, 43, 57];
    for chn = channels_num


        x_axis = -199:2:400; % in ms
        figure;
        plot(x_axis,mean_highHR(chn, :), 'color', [0 158 115]./255,  'LineWidth',1.5)
        hold on
        jbfill(x_axis, [mean_highHR(chn, :) + se_highHR(chn, :)], [mean_highHR(chn, :) - se_highHR(chn, :)], [0 158 115]./255, [0 158 115]./255, 1, 0.2)
        hold on
        plot(x_axis, mean_lowHR(chn, :), '--', 'color', [204 121 167]./255,  'LineWidth',1.5) 
        hold on
        jbfill(x_axis, [mean_lowHR(chn, :) + se_lowHR(chn, :)], [mean_lowHR(chn, :) - se_lowHR(chn, :)],[204 121 167]./255, [204 121 167]./255, 1, 0.2)
        hold on
        plot(x_axis, zeros(length(x_axis), 1), '--k',  'LineWidth', 1) 
        hold on

        % horizontal line marking the cluster showing an effect of heart
        % rate
        mask_dir = [dir_HR, filesep, groups{grp}, '\gng\'];
        load([mask_dir, 'mask_covariate.mat'])
        plot(x_axis(find(mask(chn, :) > 0) + 125), ones(length(find(mask(chn, :) == 1)), 1) * -1.25, 'k', 'LineWidth', 4.5)
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
        title([groups{grp}, ' ', name_channels] , 'FontSize', 32, 'FontWeight','normal')

    end
end


%% average across tasks
