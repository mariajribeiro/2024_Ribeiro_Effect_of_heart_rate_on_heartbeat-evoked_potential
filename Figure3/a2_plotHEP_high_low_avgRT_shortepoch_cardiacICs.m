
% plot HEP in channels FCz and CPz for participants wuth high and low
% reaction time
clear; close all
tasks = {'simpleRT', 'gng'}; group = {'young', 'older'};
load('E:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\LIMO_stats\expected_chanlocs_both.mat');

% load RT data
load('E:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\RT_median.mat');

folder_list = dir(pwd);
for grp = 1:2

    s = 0; RT = []; Yavg = [];
    for f = 1:length(folder_list)
        if contains(folder_list(f).name, 'AB')
            subj_idx = find(RT_median{grp, 1} == str2num(folder_list(f).name(3:end)));
            if ~isempty(subj_idx)
                s = s + 1;
                RT(s, 1) = mean(RT_median{grp, 2}(subj_idx, :), 2);
                load([folder_list(f).folder, filesep, folder_list(f).name, filesep, 'HEP_T_locked_data.mat'])
                Yavg(:, :, s) = mean(HEP_T_locked_data, 3);
                subjects{s} = [folder_list(f).name];
            end
        end
    end
    HEPavg_highRT{grp} = Yavg(:, :, RT > median(RT));
    HEPavg_lowRT{grp} = Yavg(:, :, RT <= median(RT));
end

%% average across tasks
groups = {'young', 'older'};
dir_RT = 'E:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\LIMO_stats\HEP_TLockedvsIBI_nobaseline_05HzFilter_cardiacICs\level2_regression_RT_avgtasks_50-350ms';

x_axis = 51:2:350; % in ms
for grp = 1:2
    mean_highRT = mean(HEPavg_highRT{grp}, 3);
    se_highRT = std(HEPavg_highRT{grp}, [], 3)/sqrt(size(HEPavg_highRT{grp}, 3));

    mean_lowRT = mean(HEPavg_lowRT{grp}, 3);
    se_lowRT = std(HEPavg_lowRT{grp}, [], 3)/sqrt(size(HEPavg_lowRT{grp}, 3));

    channels_num = 31%23;%[23 32 41];% 14%[7 16 32 34 41 50 51];%[7, 25, 43, 57];
    for chn = channels_num



        figure;
%         % grey background in time window not used for stats
%         for idx = -185:5:50
%             plot([idx idx], [-.95 1], 'color', [.9 .9 .9], 'LineWidth', 5)
%             hold on
%         end
        
        plot(x_axis,mean_highRT(chn, :), 'color', [0 158 115]./255,  'LineWidth',1.5)
        hold on
        jbfill(x_axis, [mean_highRT(chn, :) + se_highRT(chn, :)], [mean_highRT(chn, :) - se_highRT(chn, :)], [0 158 115]./255, [0 158 115]./255, 1, 0.2)
        hold on
        plot(x_axis, mean_lowRT(chn, :), '--', 'color', [204 121 167]./255,  'LineWidth',1.5) 
        hold on
        jbfill(x_axis, [mean_lowRT(chn, :) + se_lowRT(chn, :)], [mean_lowRT(chn, :) - se_lowRT(chn, :)],[204 121 167]./255, [204 121 167]./255, 1, 0.2)
        hold on
        plot(x_axis, zeros(length(x_axis), 1), '--k',  'LineWidth', 1) 
        hold on

        % horizontal line marking the cluster showing an effect of heart
        % rate
        mask_dir = [dir_RT, filesep, groups{grp}];
        if exist([mask_dir, '\mask_covariate.mat'])
            load([mask_dir, '\mask_covariate.mat'])
            plot(x_axis(find(mask(chn, :) > 0)), ones(length(find(mask(chn, :) >0)), 1) * -.8, 'k', 'LineWidth', 4.5)
            hold on
        end
        


        hold off
        box off
        ax = gca;
        axis([-inf 350 -1 1])
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
dir_RT = 'E:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\LIMO_stats\HEP_TLockedvsIBI_nobaseline_05HzFilter\level2_regression_HR';
for grp = 1:2
    mean_highRT = mean(HEP_highHR{grp, 2}, 3);
    se_highRT = std(HEP_highHR{grp, 2}, [], 3)/sqrt(size(HEP_highHR{grp, 2}, 3));

    mean_lowRT = mean(HEP_lowHR{grp, 2}, 3);
    se_lowRT = std(HEP_lowHR{grp, 2}, [], 3)/sqrt(size(HEP_lowHR{grp, 2}, 3));

    channels_num = [7 16 34 51];%[7, 25, 43, 57];
    for chn = channels_num


        figure;
        plot(x_axis,mean_highRT(chn, :), 'color', [0 158 115]./255,  'LineWidth',1.5)
        hold on
        jbfill(x_axis, [mean_highRT(chn, :) + se_highRT(chn, :)], [mean_highRT(chn, :) - se_highRT(chn, :)], [0 158 115]./255, [0 158 115]./255, 1, 0.2)
        hold on
        plot(x_axis, mean_lowRT(chn, :), '--', 'color', [204 121 167]./255,  'LineWidth',1.5) 
        hold on
        jbfill(x_axis, [mean_lowRT(chn, :) + se_lowRT(chn, :)], [mean_lowRT(chn, :) - se_lowRT(chn, :)],[204 121 167]./255, [204 121 167]./255, 1, 0.2)
        hold on
        plot(x_axis, zeros(length(x_axis), 1), '--k',  'LineWidth', 1) 
        hold on

        % horizontal line marking the cluster showing an effect of heart
        % rate
        mask_dir = [dir_RT, filesep, groups{grp}, '\gng\'];
        load([mask_dir, 'mask_covariate.mat'])
        plot(x_axis(find(mask(chn, :) > 0)), ones(length(find(mask(chn, :) == 1)), 1) * -1.25, 'k', 'LineWidth', 4.5)
        hold on
        % plot horizontal line where the groups differ - mask starts on .1
        % after R peak = .250 ms = 125 points
    %     plot(x_axis(mask(chn, :) > 0) + 250, ones(length(find(mask(chn, :) > 0)), 1) * -1.15, 's', 'MarkerSize', 4.5, 'color', 'k', ...
    %         'MarkerFaceColor', 'k')

        hold off
        box off
        ax = gca;
        axis([-inf 350 -1.5 1.3])
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
