% run second level LIMO stats
clear; close all
eeglab; limo_eeg;
% check correlation between avg ECG and avg HR - one problem could be that
% ECG electrodes position was not well controled across participants
% addpath('C:\eeglab2022.0\plugins\limo_tools-master_Nov2022') 

older = [11 12 14 20 21 22 23 32 37 38 41 43 47 48 49 52 55 57 58 63 64 65 67 69 7 70 71 75 8 83 86];
young = [13 15 16 25 26 28 31 33 34 36 4 42 44 45 46 50 51 53 54 56 59 6 62 66 68 72 74 76 78 80 82 84 85 9];
tasks = {'simpleRT', 'gng'}; group = {'young', 'older'};
folder_list = dir(pwd);
load('D:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\LIMO_stats\expected_chanlocs_both.mat');

% load HR data calculated in D:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\hep_effect_of_heart_rate.m
load('D:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\HR_avg_std_min_max_range.mat');


%% limo_random_robust(4,y,X,parameter number,LIMO)
%                    4 = regression analysis
%                    y = data (dim channels, time or freq, subjects)
%                      = data (dim channels, freq, time, subjects)
%                    X = continuous regressor(s)
%                    parameter number = describe which parameter is currently analysed (e.g. 1 - use for maming only)

for grp = 1:2

        level2_regression_HR_dir = [pwd, '\level2_regression_ECG_HR_avgtasks_50-350ms\', group{grp}];
        mkdir(level2_regression_HR_dir);
        cd(level2_regression_HR_dir);
        
        s = 0; Y = []; HR = [];
        for f = 1:length(folder_list)
            if contains(folder_list(f).name, 'AB')
                % check if participant has HR data
                subj_idx = find(HR_avg_std_min_max_range{grp, 1} == str2num(folder_list(f).name(3:end)));
                if ~isempty(subj_idx)
                    s = s + 1;
                    HR(s, 1) = mean([HR_avg_std_min_max_range{grp, 2}(subj_idx, 2, 1), ...
                       HR_avg_std_min_max_range{grp, 2}(subj_idx, 3, 1)]);
                    load([folder_list(f).folder, filesep, folder_list(f).name, filesep, 'ECG_T_locked_data.mat'])
                    Y(1, :, s) = mean(ECG_T_locked_data(1, 1:150, :), 3); %channels, time, subjects)
                    subjects{s} = [folder_list(f).name];
                end
            end
        end

        save Y Y
        save HR HR

        % create LIMO variable
        LIMO.Level                    = 2;

        LIMO.Type = 'Channels';
        LIMO.Analysis = 'Time';

        LIMO.dir                      = level2_regression_HR_dir;

        LIMO.data.data_dir            = level2_regression_HR_dir;
        % LIMO.data.data                = 

        LIMO.data.start               = 51;
        LIMO.data.end                 = 349;
        LIMO.data.trim1               = 1; 
        LIMO.data.trim2               = size(Y, 2);
        LIMO.data.timevect            = 51:2:349;
        LIMO.data.sampling_rate       = 500;
        LIMO.data.neighbouring_matrix = channeighbstructmat(1, 2);
        LIMO.data.chanlocs = expected_chanlocs(1);

        % LIMO.design.fullfactorial     = 0;
        % LIMO.design.zscore            = 0;
    %     LIMO.design.method            = 'Trimmed mean'; 
        LIMO.design.type_of_analysis  = 'Mass-univariate';
        LIMO.design.bootstrap         = 1000;
        LIMO.design.tfce              = 0;

        save LIMO LIMO

        % regression

        % limo_random_robust(4,y,X,parameter number,LIMO)
        LIMOpath = limo_random_robust(4, Y, HR, 1,LIMO);


        %% view results  - clustering
        % addpath('C:\eeglab2022.0\plugins\limo_tools-master_Nov2022') 
        % level2_regression_pupil_dir = [pwd, '\level2_RM_ANOVA_IBIvsHEP'];
        % cd(level2_regression_pupil_dir);
        % https://github.com/LIMO-EEG-Toolbox/limo_tools/wiki/Workspace_variables
        mask = [];
        % limo_display_results(Type,FileName,PathName,p,MCC,LIMO,flag,options)
        limo_display_results(1,'Covariate_effect_1.mat', pwd,0.05, 2,...
            fullfile(pwd,'LIMO.mat'), 0);

        if ~isempty(mask)
            save mask_covariate mask
            save stat_covariate stat_values %F/t values
            save p_values_covariate p_values


            number_clusters = max(unique(mask(:)));
            p_value_cluster = p_values(mask >0 );
            max_p_value = max(p_values(mask ~= 0));
            %     what is the range of significant F/t values? 
            stats_F = [min(stat_values(mask ~= 0)) max(stat_values(mask ~= 0))];
        end
        %%
        cd ..
        cd ..

end


%% determine electrodes belonging to cluster and time window
tasks = {'simpleRT', 'gng'}; group = {'young', 'older'};
chn_clust = cell(2);
tim_clust = cell(2);
chn_names_clust = cell(2);
load('D:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\LIMO_stats\expected_chanlocs_both.mat');

for grp = 1:2
    level2_regression_HR_dir = [pwd, '\level2_regression_ECG_HR_avgtasks_50-350ms\', group{grp}];
        cd(level2_regression_HR_dir);
        if isfile('mask_covariate.mat')
            load('mask_covariate.mat')

            [chn, time] = find(mask > 0);

            chn_clust{grp} = unique(chn);
            
            tim_clust{grp} = unique(time);
            ch=0;
            for c = chn_clust{grp}'
                ch = ch + 1;
                chn_names_clust{grp}{ch} = expected_chanlocs(c).labels;
            end
        end
        
        cd ..
        cd ..
        
end

%% plot avg HEP within sig cluster vs heart rate
% older group
level2_regression_HR_dir = [pwd, '\level2_regression_ECG_HR_avgtasks_50-350ms\older'];

load([level2_regression_HR_dir, filesep, 'Y']); 
load([level2_regression_HR_dir, filesep, 'mask_covariate']); 
load([level2_regression_HR_dir, filesep, 'HR']);

[row,col] = find(mask > 0);
ECG_sig = [];
for s = 1:size(Y, 3)
    for idx = 1:length(row)
        ECG_sig(s, idx) = Y(row(idx), col(idx), s);
    end
end

ECG_avg = mean(ECG_sig, 2)


[R_older,P_older] = corrcoef(ECG_avg, HR)


figure;
plot(HR, ECG_avg,'d', 'color', 'k',  'MarkerSize',12, 'LineWidth', 1.5); 
if P_older(1, 2)<.05
    lsline;
end
ax = gca;   box off
ax.LineWidth = 2.5;
ax.FontSize = 28;
ax.FontName = 'Arial';
ylabel({'ECG amplitude', '(\muV)'}, 'FontSize', 32, 'FontWeight','normal');
xlabel('Heart rate (bpm)', 'FontSize', 32, 'FontWeight','normal');
title('Older', 'FontSize', 32, 'FontWeight','normal');
% axis([250 500 -inf inf])
ylim padded
xlim padded


%% young group
level2_regression_HR_dir = [pwd, '\level2_regression_ECG_HR_avgtasks_50-350ms\young'];

load([level2_regression_HR_dir, filesep, 'Y']); 
% load([level2_regression_HR_dir, filesep, 'mask_covariate']); 
load([level2_regression_HR_dir, filesep, 'HR']);

ECG_sig = [];
for s = 1:size(Y, 3)
    for idx = 1:length(row)
        ECG_sig(s, idx) = Y(row(idx), col(idx), s);
    end
end

ECG_avg = mean(ECG_sig, 2)



[R_young,P_young] = corrcoef(ECG_avg, HR)

ECG = ECG_avg(ECG_avg > -100);
HR_clean = HR(ECG_avg > -100);
[R_yng,P_yng] = corrcoef(ECG, HR_clean)

figure;
plot(HR, ECG_avg,'d', 'color', 'k',  'MarkerSize',12, 'LineWidth', 1.5); 
if P_older(1, 2)<.05
    lsline;
end
ax = gca;  box off
ax.LineWidth = 2.5;
ax.FontSize = 28;
ax.FontName = 'Arial';
ylabel({'ECG amplitude', '(\muV)'}, 'FontSize', 32, 'FontWeight','normal');
xlabel('Heart rate (bpm)', 'FontSize', 32, 'FontWeight','normal');
title('Young', 'FontSize', 32, 'FontWeight','normal');
% axis([250 500 -inf inf])
ylim padded
xlim padded


figure;
plot(HR_clean, ECG,'d', 'color', 'k',  'MarkerSize',12, 'LineWidth', 1.5); 
if P_older(1, 2)<.05
    lsline;
end
ax = gca;   box off
ax.LineWidth = 2.5;
ax.FontSize = 28;
ax.FontName = 'Arial';
ylabel({'ECG amplitude', '(\muV)'}, 'FontSize', 32, 'FontWeight','normal');
xlabel('Heart rate (bpm)', 'FontSize', 32, 'FontWeight','normal');
title('Young', 'FontSize', 32, 'FontWeight','normal');
% axis([250 500 -inf inf])
ylim padded
xlim padded

%% highest F values
level2_regression_HR_dir = [pwd, '\level2_regression_HR_avgtasks_50-350ms\older'];
load([level2_regression_HR_dir, filesep, 'Covariate_effect_1']); 
F_mean = mean(Covariate_effect(:, :, 1), 2)

[B,index] = sortrows(F_mean)
c=0; clear chn_F
for ch = index'
    c = c+1;
    chn_F{c} = expected_chanlocs(ch).labels;
end

%% topo plot
group = {'young', 'older'};
for grp = 1:2
    level2_regression_HR_dir = [pwd, '\level2_regression_HR_avgtasks_50-350ms\', group{grp}];

    % plot topographic plot of F-values with significant channels highlighted
    % channels where the T-HEP is significantly different from zero
    mask = [];
    if exist([level2_regression_HR_dir, filesep 'mask_covariate.mat'])
        load([level2_regression_HR_dir, filesep 'mask_covariate']);
    end
    
    load([level2_regression_HR_dir, '/Covariate_effect_1']);

    F_value = squeeze(mean(Covariate_effect(:, 100:end, 1), 2));
    
    sig_chan_number = find(sum(mask, 2) > 0);

    title_txt = {'Young' 'Older'};

    figure;
    if~isempty(sig_chan_number)
        topoplot(ones(59, 1)*.025, expected_chanlocs, 'electrodes', 'off', 'plotchans', sig_chan_number, 'style', 'blank',...
            'plotdisk', 'on',  'hcolor'  , 'none') ; hold on
    end
    topoplot(F_value, expected_chanlocs, 'electrodes', 'off'); 
    caxis([-0 10]); %c.Axis.FontSize = 16;
    cH = colorbar; set(cH,'FontSize',30);
    % colorbar('Ticks',[0, 1, 2, 3, 4, 5], 'FontSize', 18, 'FontWeight','normal');
    colormap(crameri('imola')); % needs colour maps from http://www.fabiocrameri.ch/colourmaps.php
    % title('Young', 'FontSize', 30, 'FontWeight','normal')
    % set(get(gca,'title'),'Position',[0,-.65, 0])
    %         text(5, 0.4, title_txt{task, grp})
    
end