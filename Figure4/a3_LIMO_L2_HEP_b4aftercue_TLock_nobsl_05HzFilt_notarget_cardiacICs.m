% LIMO analyses comparing before and after cue in simple RT and gng task in
% both age groups - only Level 2 as the correction for baseline shift does
% not work on single trial. for each participant we only have average
% across trials - input that for Level 2
clear; close all
% smooth data to subtract from HEP
files_dir = 'E:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\TLockedNoBaseline\NoBaseline05HzFilt_NoTarget_cardiacICs';
tasks = {'passive', 'simpleRT', 'gonogo'};
load([files_dir, filesep, 'HEP_T_locked_aftercue_random']); 
load([files_dir, filesep, 'HEP_T_locked_b4cue_random']); 
load([files_dir, filesep, 'HEP_T_locked_random']); 

% load HEPs
load([files_dir, filesep, 'HEP_T_locked_aftercue']); 
load([files_dir, filesep, 'HEP_T_locked_b4cue']); 
load([files_dir, filesep, 'HEP_T_locked_all']);

HEP_T_aftercue_detrended = HEP_T_locked_aftercue;
HEP_T_b4cue_detrended = HEP_T_locked_b4cue;


smth_HEP_T_locked_aftercue_random = HEP_T_locked_aftercue;
smth_HEP_T_locked_b4cue_random = HEP_T_locked_b4cue;
% smth_HEP_T_locked_random = HEP_T_locked_random;

smth_window_length = 150;

for grp = 1:2
    for p = 1:size(HEP_T_locked_b4cue{grp, 2}, 1)
        for t = 1:3
            subject = HEP_T_locked_aftercue_random{grp, 1}(p); % should b the same for HEP_T_locked_b4cue_random
            % smooth ERP obtained with random events after cue
            data_tmp = squeeze(HEP_T_locked_aftercue_random{grp, 2}(p, t, :, :)); % channels x time
            data_tmp_smooth = smoothdata(data_tmp, 2, 'gaussian', smth_window_length);
            smth_HEP_aftercue_random{grp, 2}(p, t, :, :) = data_tmp_smooth;
            
            % subtract smoothed ERP from after cue HEP to correct for
            % baseline shift
            data = [];
            for sbj = 1:size(HEP_T_locked_aftercue{grp, 2}, 1)
                if HEP_T_locked_aftercue{grp, 1}(sbj) == subject
                    data = squeeze(HEP_T_locked_aftercue{grp, 2}(sbj, t, :, :)); % channels x time
                end
            end
            HEP_T_aftercue_detrended{grp, 2}(p, t, :, :) = data - data_tmp_smooth;
            
            % smooth HEP-random before cue
            data_tmp = squeeze(HEP_T_locked_b4cue_random{grp, 2}(p, t, :, :)); % channels x time
            data_tmp_smooth = smoothdata(data_tmp, 2, 'gaussian', smth_window_length);
            smth_HEP_b4cue_random{grp, 2}(p, t, :, :) = data_tmp_smooth;
            
            % correct for baseline shift HEP before cue
            data = [];
            for sbj = 1:size(HEP_T_locked_b4cue{grp, 2}, 1)
                if HEP_T_locked_b4cue{grp, 1}(sbj) == subject
                    data = squeeze(HEP_T_locked_b4cue{grp, 2}(sbj, t, :, :)); % channels x time
                end
            end
            HEP_T_b4cue_detrended{grp, 2}(p, t, :, :) = data - data_tmp_smooth;
            
            data_tmp = squeeze(HEP_T_locked_random{grp, 2}(p, t, :, :)); % channels x time
            data_tmp_smooth = smoothdata(data_tmp, 2, 'gaussian', smth_window_length);
            smth_HEP_random{grp, 2}(p, t, :, :) = data_tmp_smooth;
        end
    end
end

eeglab; limo_eeg;
load('E:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\LIMO_stats\expected_chanlocs_both.mat');
%%
% repeated measures ANOVA within-subject factors before cue x after cue; simple RT x go/no-go
% limo_random_robust(6,y,gp,factor_levels,LIMO,'go',option)
%                    6 = 

%                    y = data (dim channels, time or freq, subjects, measures)
%                      = data (dim channels, freq, time, subjects, measures)
%                      = the name of the Yr file
%                    gp = a vector defining gps
%                    factor_levels = a vector specifying the levels of each repeated measure factor
%                    LIMO the basic structure with data, design and channel info
%                         or the full name of the LIMO file
%                    'go' is optional and prompt or not the pseudo-design
%                         options are 'yes' (prompt, default),
%                         or 'no' (no prompt, usuful for scripting)

%  y = data (dim channels, time or freq, subjects, measures)
% I think the measures should be task x timing = simple RT before cue,
% simple RT after cue, go/no-go before cue, go/no-go after cue


% HEP_T_aftercue_detrended{grp, 2}(p, t, :, :) =  participant x task x
% channel x time

y = NaN(59, length(126:300), length(HEP_T_b4cue_detrended{1, 1})+length(HEP_T_b4cue_detrended{2, 1}), 4);
gp = [];

y_tmp{1} = NaN(59, length(176:325), length(HEP_T_b4cue_detrended{1, 1}), 4);
y_tmp{2} = NaN(59, length(176:325), length(HEP_T_b4cue_detrended{2, 1}), 4);


for grp = 1:2
    for s = 1:size(HEP_T_b4cue_detrended{grp, 2}, 1)
        % simple RT before cue
        y_tmp{grp}(:, :, s, 1) = HEP_T_b4cue_detrended{grp, 2}(s, 2, 1:59, 176:end); % analyse only from 51ms onwards
        % simple RT after cue
        y_tmp{grp}(:, :, s, 2) = HEP_T_aftercue_detrended{grp, 2}(s, 2, 1:59, 176:end);
        % gng before cue
        y_tmp{grp}(:, :, s, 3) = HEP_T_b4cue_detrended{grp, 2}(s, 3, 1:59, 176:end);
        % gng after cue
        y_tmp{grp}(:, :, s, 4) = HEP_T_aftercue_detrended{grp, 2}(s, 3, 1:59, 176:end);
        
        gp = [gp; grp];
    end
end

% figure; plot(squeeze(mean(y_tmp{1}(37, :, :, 4), 3)))

y = cat(3, y_tmp{1}, y_tmp{2});


%% eeglab; limo_eeg;
load('E:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\LIMO_stats\expected_chanlocs_both.mat');

level2_RM_ANOVA_dir = [pwd filesep 'level2_RM_ANOVA_HEP_TLocked_b4_after_cue'];
mkdir(level2_RM_ANOVA_dir)
cd(level2_RM_ANOVA_dir)

eeglab;
% load('E:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\LIMO_stats\expected_chanlocs_both.mat');

% create LIMO variable
LIMO.Level                    = 2;

LIMO.Type = 'Channels';
LIMO.Analysis = 'Time';

LIMO.dir                      = level2_RM_ANOVA_dir;

LIMO.data.data_dir            = level2_RM_ANOVA_dir;
% LIMO.data.data                = 
LIMO.data.start               = 51;
LIMO.data.end                 = 349;
LIMO.data.trim1               = 1; 
LIMO.data.trim2               = size(y, 2);
LIMO.data.timevect            = 51:2:349;
LIMO.data.sampling_rate       = 500;
LIMO.data.neighbouring_matrix = channeighbstructmat;
LIMO.data.chanlocs = expected_chanlocs;

% LIMO.design.fullfactorial     = 0;
% LIMO.design.zscore            = 0;
LIMO.design.method            = 'Trimmed mean'; 
LIMO.design.type_of_analysis  = 'Mass-univariate';
LIMO.design.bootstrap         = 1000; 
LIMO.design.tfce              = 0;

save LIMO LIMO

factor_levels = [2, 2];
%
limo_random_robust(6,y,gp,factor_levels,LIMO)

%% view results  - clustering
addpath('C:\eeglab2022.0\plugins\limo_tools-master_Nov2022') 
% level2_RM_ANOVA_dir = [pwd, '\level2_RM_ANOVA_IBIvsHEP'];
% cd(level2_RM_ANOVA_dir);
% https://github.com/LIMO-EEG-Toolbox/limo_tools/wiki/Workspace_variables
mask = [];
% limo_display_results(Type,FileName,PathName,p,MCC,LIMO,flag,options)
limo_display_results(1,'Rep_ANOVA_Main_effect_1.mat',pwd,0.1, 2,...
    fullfile(pwd,'LIMO.mat'), 0);

if ~isempty(mask)
    save mask_main_effect mask
    save stat_values_main_effect stat_values %F/t values
    save p_values_main_effect p_values


    number_clusters = max(unique(mask(:)));
    p_value_cluster = p_values(mask == 1);
    max_p_value = max(p_values(mask ~= 0));
    %     what is the range of significant F/t values? 
    stats_F = [min(stat_values(mask ~= 0)) max(stat_values(mask ~= 0))];
end

%%
mask = [];
limo_display_results(1,'Rep_ANOVA_Interaction_gp_Factor_1.mat',pwd,0.05, 2,...
    fullfile(pwd,'LIMO.mat'), 0);

%%
mask = [];
limo_display_results(1,'Rep_ANOVA_Interaction_gp_Factor_2.mat',pwd,0.05, 2,...
    fullfile(pwd,'LIMO.mat'), 0);

%%
mask = [];
limo_display_results(1,'Rep_ANOVA_Interaction_gp_Factors_12.mat',pwd,0.05, 2,...
    fullfile(pwd,'LIMO.mat'), 0);

%%
mask = [];
limo_display_results(1,'Rep_ANOVA_Gp_effect.mat',pwd,0.1, 2,...
    fullfile(pwd,'LIMO.mat'), 0);
% saveas(gcf, 'One_sample_timecourse.fig'); close(gcf)
if ~isempty(mask)
    save mask_group_effect mask
    save stat_values_group_effect stat_values %F/t values
    save p_values_group_effect p_values
    
    number_clusters = max(unique(mask(:)));
    p_value_cluster = p_values(mask == 1);
    max_p_value = max(p_values(mask ~= 0));
    %     what is the range of significant F/t values? 
    stats_F = [min(stat_values(mask ~= 0)) max(stat_values(mask ~= 0))];
    
    summary_stats = limo_get_summary('Rep_ANOVA_Gp_effect.mat',mask);
end

%%
mask = [];
limo_display_results(1,'Rep_ANOVA_Main_effect_2.mat',pwd,0.1, 2,...
    fullfile(pwd,'LIMO.mat'), 0);
% saveas(gcf, 'One_sample_timecourse.fig'); close(gcf)
if ~isempty(mask)
    save mask_main_effect_2 mask
    save stat_values_main_effect_2 stat_values %F/t values
    save p_values_main_effect_2 p_values
    
    number_clusters = max(unique(mask(:)))
    p_value_cluster = p_values(mask == 1);
    max_p_value = max(p_values(mask ~= 0))
    %     what is the range of significant F/t values? 
    stats_F = [min(stat_values(mask ~= 0)) max(stat_values(mask ~= 0))];
    
    summary_stats = limo_get_summary('Rep_ANOVA_Main_effect_2.mat',mask);
end

%%
% use only paired t-test - all participants together
cd('E:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\LIMO_stats\HEP_b4_after_cue_TLocked_nobaslin_05HzFilter_notarget');

level2_onesample_dir = [pwd, '\level2_b4cue_vs_aftercue'];
tasks = {'passive', 'simpleRT', 'gonogo'};
% separate analyses per task - simple RT and gng
for t = 2:3
    
    % input for stats T-HEP in the interval 50-400 ms after T peaks
    % data for LIMO paired t-test b4cue vs after cue - all participants together
    % y1 = b4cue - HEP_T_locked_b4cue_detrended{grp, 2}(sbj, t, :, :) - data (dim channels, time, subjects)
    % y2 = after cue - HEP_T_locked_aftercue_detrended{grp, 2}(sbj, t, :, :)
    y1_yng = []; y2_yng = [];
    for s = 1:size(HEP_T_b4cue_detrended{1, 2}, 1)
        y1_yng(:, :, s) = squeeze(HEP_T_b4cue_detrended{1, 2}(s, t, 1:59, 126:300));
        y2_yng(:, :, s) = squeeze(HEP_T_aftercue_detrended{1, 2}(s, t, 1:59, 126:300));
    end
    
    y1_old = []; y2_old = [];
    for s = 1:size(HEP_T_b4cue_detrended{2, 2}, 1)
        y1_old(:, :, s) = squeeze(HEP_T_b4cue_detrended{2, 2}(s, t, 1:59, 126:300));
        y2_old(:, :, s) = squeeze(HEP_T_aftercue_detrended{2, 2}(s, t, 1:59, 126:300));
    end
    
    y1 = cat(3, y1_yng, y1_old);
    y2 = cat(3, y2_yng, y2_old);

    % limo_random_robust(3,y1,y2,parameter number,LIMO)
    %                    3 = paired t-test
    %                    y1 = data (dim channels, time or freq, subjects)
    %                       = data (dim channels, freq, time, subjects)
    %                       = the name of the Y1r file
    %                    y2 = data (dim channels, time or freq, subjects)
    %                       = data (dim channels, freq, time, subjects)
    %                       = the name of the Y2r file
    %                    parameter number = describe which parameter is currently analysed (e.g. 1 - use for maming only)
    
    LIMO_directory = [level2_onesample_dir filesep tasks{t}];
    mkdir(LIMO_directory);
    cd(LIMO_directory);
    
    % create LIMO variable
    LIMO.Level                    = 2;

    LIMO.Type = 'Channels';
    LIMO.Analysis = 'Time';

    LIMO.dir                      = LIMO_directory;

    LIMO.data.data_dir            = LIMO_directory;
    % LIMO.data.data                = 

    LIMO.data.start               = 51;
    LIMO.data.end                 = 400;
    LIMO.data.trim1               = 1; 
    LIMO.data.trim2               = size(y1, 2);
    LIMO.data.timevect            = 51:2:400;
    LIMO.data.sampling_rate       = 500;
    LIMO.data.neighbouring_matrix = channeighbstructmat;
    LIMO.data.chanlocs = expected_chanlocs;

    % LIMO.design.fullfactorial     = 0;
    % LIMO.design.zscore            = 0;
    LIMO.design.method            = 'Trimmed mean'; 
    LIMO.design.type_of_analysis  = 'Mass-univariate';
    LIMO.design.bootstrap         = 1000;  
    LIMO.design.tfce              = 0;

    save LIMO LIMO

%     limo_random_robust(3,y1,y2,parameter number,LIMO)
    limo_random_robust(3, y1, y2, 1, LIMO)

    %% view results  - clustering
    % https://github.com/LIMO-EEG-Toolbox/limo_tools/wiki/Workspace_variables
    mask = [];
    % limo_display_results(Type,FileName,PathName,p,MCC,LIMO,flag,options)
    limo_display_results(1,'paired_samples_ttest_parameter_1', pwd, 0.05, 2,...
        fullfile(pwd,'LIMO.mat'), 0);
    % saveas(gcf, 'One_sample_timecourse.fig'); close(gcf)
    if ~isempty(mask)
        save mask mask
        save stat_values stat_values %F/t values
        save p_values p_values
        
        [row,col] = find(mask == 1);
        max_time = max(col)*2+100;

        number_clusters = max(unique(mask(:)));
        p_value_cluster = p_values(mask == 1);
        max_p_value = max(p_values(mask ~= 0));
        %     what is the range of significant F/t values? 
        stats_F = [min(stat_values(mask ~= 0)) max(stat_values(mask ~= 0))];
    end
%     % replot using 
%     limo_display_image(LIMO,stat_values ,mask,'go/no-go HR effect',0) 
end

%% b4cue vs after cue
% use only paired t-test - all participants together - all tasks together
cd('E:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\LIMO_stats\HEP_b4_after_cue_TLocked_nobaslin_05HzFilter_notarget');

level2_onesample_dir = [pwd, '\level2_b4cue_vs_aftercue'];
tasks = {'passive', 'simpleRT', 'gonogo'};
% separate analyses per task - simple RT and gng

% input for stats T-HEP in the interval 50-400 ms after T peaks
% data for LIMO paired t-test b4cue vs after cue - all participants together
% y1 = b4cue - HEP_T_locked_b4cue_detrended{grp, 2}(sbj, t, :, :) - data (dim channels, time, subjects)
% y2 = after cue - HEP_T_locked_aftercue_detrended{grp, 2}(sbj, t, :, :)
y1_yng = []; y2_yng = [];
for s = 1:size(HEP_T_b4cue_detrended{1, 2}, 1)
    y1_yng(:, :, s) = squeeze(mean(HEP_T_b4cue_detrended{1, 2}(s, 2:3, 1:59, 176:end), 2));
    y2_yng(:, :, s) = squeeze(mean(HEP_T_aftercue_detrended{1, 2}(s, 2:3, 1:59, 176:end), 2));
end

y1_old = []; y2_old = [];
for s = 1:size(HEP_T_b4cue_detrended{2, 2}, 1)
    y1_old(:, :, s) = squeeze(mean(HEP_T_b4cue_detrended{2, 2}(s, 2:3, 1:59, 176:end), 2));
    y2_old(:, :, s) = squeeze(mean(HEP_T_aftercue_detrended{2, 2}(s, 2:3, 1:59, 176:end), 2));
end

y1 = cat(3, y1_yng, y1_old);
y2 = cat(3, y2_yng, y2_old);

% limo_random_robust(3,y1,y2,parameter number,LIMO)
%                    3 = paired t-test
%                    y1 = data (dim channels, time or freq, subjects)
%                       = data (dim channels, freq, time, subjects)
%                       = the name of the Y1r file
%                    y2 = data (dim channels, time or freq, subjects)
%                       = data (dim channels, freq, time, subjects)
%                       = the name of the Y2r file
%                    parameter number = describe which parameter is currently analysed (e.g. 1 - use for maming only)

LIMO_directory = [level2_onesample_dir filesep 'avg_tasks'];
mkdir(LIMO_directory);
cd(LIMO_directory);

% create LIMO variable
LIMO.Level                    = 2;

LIMO.Type = 'Channels';
LIMO.Analysis = 'Time';

LIMO.dir                      = LIMO_directory;

LIMO.data.data_dir            = LIMO_directory;
% LIMO.data.data                = 

LIMO.data.start               = 51;
LIMO.data.end                 = 400;
LIMO.data.trim1               = 1; 
LIMO.data.trim2               = size(y1, 2);
LIMO.data.timevect            = 51:2:400;
LIMO.data.sampling_rate       = 500;
LIMO.data.neighbouring_matrix = channeighbstructmat;
LIMO.data.chanlocs = expected_chanlocs;

% LIMO.design.fullfactorial     = 0;
% LIMO.design.zscore            = 0;
LIMO.design.method            = 'Trimmed mean'; 
LIMO.design.type_of_analysis  = 'Mass-univariate';
LIMO.design.bootstrap         = 1000;  
LIMO.design.tfce              = 0;

save LIMO LIMO

%     limo_random_robust(3,y1,y2,parameter number,LIMO)
limo_random_robust(3, y1, y2, 1, LIMO)

% view results  - clustering
% https://github.com/LIMO-EEG-Toolbox/limo_tools/wiki/Workspace_variables
mask = [];
% limo_display_results(Type,FileName,PathName,p,MCC,LIMO,flag,options)
limo_display_results(1,'paired_samples_ttest_parameter_1', pwd, 0.05, 2,...
    fullfile(pwd,'LIMO.mat'), 0);
% saveas(gcf, 'One_sample_timecourse.fig'); close(gcf)
if ~isempty(mask)
    save mask mask
    save stat_values stat_values %F/t values
    save p_values p_values

    [row,col] = find(mask == 1);
    max_time = max(col)*2+100;

    number_clusters = max(unique(mask(:)));
    p_value_cluster = p_values(mask == 1);
    max_p_value = max(p_values(mask ~= 0));
    %     what is the range of significant F/t values? 
    stats_F = [min(stat_values(mask ~= 0)) max(stat_values(mask ~= 0))];
end
%     % replot using 
%     limo_display_image(LIMO,stat_values ,mask,'go/no-go HR effect',0) 


%%
% comparison across tasks
% use only paired t-test - all participants together - b4 and after cue
% together
cd('E:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\LIMO_stats\HEP_b4_after_cue_TLocked_nobaslin_05HzFilter_notarget');

level2_onesample_dir = [pwd, '\level2_simpleRT_vs_gng'];
tasks = {'passive', 'simpleRT', 'gonogo'};
% separate analyses per task - simple RT and gng

    
% input for stas T-HEP in the interval 50-400 ms after T peaks
% data for LIMO paired t-test b4cue vs after cue - all participants together
% y1 = b4cue - HEP_T_locked_b4cue_detrended{grp, 2}(sbj, t, :, :) - data (dim channels, time, subjects)
% y2 = after cue - HEP_T_locked_aftercue_detrended{grp, 2}(sbj, t, :, :)
y1_yng = []; y2_yng = [];

% average HEP per participant average b4cue and after cue
HEP_avg = cell(2, 1);
for grp = 1:2
    for s = 1:size(HEP_T_b4cue_detrended{grp, 2}, 1)
        for t = 1:3
            for tm = 1:size(HEP_T_b4cue_detrended{1, 2}, 4)
                for chn = 1:59
                    HEP_avg{grp}(s, t, chn, tm) = mean([HEP_T_b4cue_detrended{grp, 2}(s, t, chn, tm),...
                        HEP_T_aftercue_detrended{grp, 2}(s, t, chn, tm)]);
                end
            end
        end
    end
end


y1_yng_simpleRT = []; y1_yng_gng = []; 
for s = 1:size(HEP_avg{1}, 1)
    y1_yng_simpleRT(:, :, s) = squeeze(HEP_avg{1}(s, 2, 1:59, 126:300));
    y1_yng_gng(:, :, s) = squeeze(HEP_avg{1}(s, 3, 1:59, 126:300));
end
y1_old_simpleRT = [];  y1_old_gng = [];
for s = 1:size(HEP_avg{2}, 1)
    y1_old_simpleRT(:, :, s) = squeeze(HEP_avg{2}(s, 2, 1:59, 126:300));
    y1_old_gng(:, :, s) = squeeze(HEP_avg{2}(s, 3, 1:59, 126:300));
end

y_simpleRT = cat(3, y1_yng_simpleRT, y1_old_simpleRT);
y_gng = cat(3, y1_yng_gng, y1_old_gng);

% limo_random_robust(3,y1,y2,parameter number,LIMO)
%                    3 = paired t-test
%                    y1 = data (dim channels, time or freq, subjects)
%                       = data (dim channels, freq, time, subjects)
%                       = the name of the Y1r file
%                    y2 = data (dim channels, time or freq, subjects)
%                       = data (dim channels, freq, time, subjects)
%                       = the name of the Y2r file
%                    parameter number = describe which parameter is currently analysed (e.g. 1 - use for maming only)

LIMO_directory = [level2_onesample_dir filesep tasks{t}];
mkdir(LIMO_directory);
cd(LIMO_directory);

% create LIMO variable
LIMO.Level                    = 2;

LIMO.Type = 'Channels';
LIMO.Analysis = 'Time';

LIMO.dir                      = LIMO_directory;

LIMO.data.data_dir            = LIMO_directory;
% LIMO.data.data                = 

LIMO.data.start               = 51;
LIMO.data.end                 = 400;
LIMO.data.trim1               = 1; 
LIMO.data.trim2               = size(y_simpleRT, 2);
LIMO.data.timevect            = 51:2:400;
LIMO.data.sampling_rate       = 500;
LIMO.data.neighbouring_matrix = channeighbstructmat;
LIMO.data.chanlocs = expected_chanlocs;

% LIMO.design.fullfactorial     = 0;
% LIMO.design.zscore            = 0;
LIMO.design.method            = 'Trimmed mean'; 
LIMO.design.type_of_analysis  = 'Mass-univariate';
LIMO.design.bootstrap         = 1000;  
LIMO.design.tfce              = 0;

save LIMO LIMO

%     limo_random_robust(3,y1,y2,parameter number,LIMO)
limo_random_robust(3, y_simpleRT, y_gng, 1, LIMO)

% view results  - clustering
% https://github.com/LIMO-EEG-Toolbox/limo_tools/wiki/Workspace_variables
mask = [];
% limo_display_results(Type,FileName,PathName,p,MCC,LIMO,flag,options)
limo_display_results(1,'paired_samples_ttest_parameter_1', pwd, 0.05, 2,...
    fullfile(pwd,'LIMO.mat'), 0);
% saveas(gcf, 'One_sample_timecourse.fig'); close(gcf)
if ~isempty(mask)
    save mask mask
    save stat_values stat_values %F/t values
    save p_values p_values

    [row,col] = find(mask == 1);
    max_time = max(col)*2+100;

    number_clusters = max(unique(mask(:)));
    p_value_cluster = p_values(mask == 1);
    max_p_value = max(p_values(mask ~= 0));
    %     what is the range of significant F/t values? 
    stats_F = [min(stat_values(mask ~= 0)) max(stat_values(mask ~= 0))];
end
%     % replot using 
%     limo_display_image(LIMO,stat_values ,mask,'go/no-go HR effect',0) 


%% plot b4cue and after cue   - both groups together - simple RT vs gng together
% tasks together

b4cue_data = []; aftercue_data =[];
b4cue_data = squeeze(mean(y(:, :, :, [1, 3]), 4));
aftercue_data = squeeze(mean(y(:, :, :, [2, 4]), 4));

% HEP b4cue
b4cue_mean = squeeze(mean(b4cue_data, 3, 'omitnan'));
b4cue_se = squeeze(std(b4cue_data, [], 3, 'omitnan'))/sqrt(size(b4cue_data, 3));

% HEP aftercue
aftercue_mean = squeeze(mean(aftercue_data, 3, 'omitnan'));
aftercue_se = squeeze(std(aftercue_data, [], 3, 'omitnan'))/sqrt(size(aftercue_data, 3));

chan_dir = 'E:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study';
channels = load([chan_dir filesep 'chanlocs.mat']);
channels_num = [7, 16, 25, 34, 51];

for p = channels_num

    x_axis=51:2:400; % in mseconds
    figure;
    plot(x_axis, b4cue_mean(p, :), 'color', [0 0.4470 0.7410],  'LineWidth',1.5) %plot mean gng (light blue)
    hold on
    jbfill(x_axis, [b4cue_mean(p, :) + b4cue_se(p, :)], [b4cue_mean(p, :) - b4cue_se(p, :)], [0 0.4470 0.7410],[0 0.4470 0.7410], 1, 0.2)
    hold on
    plot(x_axis, aftercue_mean(p, :), 'color', 'k',  'LineWidth',1.5) %plot mean gng (light blue)
    hold on
    jbfill(x_axis, [aftercue_mean(p, :) + aftercue_se(p, :)], [aftercue_mean(p, :) - aftercue_se(p, :)], 'k', 'k', 1, 0.2)
    hold on
    plot(x_axis, zeros(length(x_axis), 1), ':k')
    hold off
    box off
    ax = gca;
    if t == 1
        axis([-inf inf -1.5 1.5])
    else
        axis([-inf inf -2.5 1.5])
    end
    ax.LineWidth = 2.5;
    ax.FontSize = 28;
    ax.FontName = 'Arial';
    xlabel('Time (ms)', 'FontSize', 32, 'FontWeight','normal')
    ylabel('Amplitude (\muV)', 'FontSize', 32, 'FontWeight','normal')
    name_channels =  channels.chanlocs(p).labels;
    title(name_channels, 'FontSize', 32, 'FontWeight','normal')
%     legend(' before cue', 'b4', 'after cue', 'after')
end

%% plot simple RT vs gng - both groups together - b4cue and after cue together
% tasks together

simpleRT_data = []; gng_data =[];
simpleRT_data = squeeze(mean(y(:, :, :, 1:2), 4));
gng_data = squeeze(mean(y(:, :, :, 3:4), 4));

% HEP simpleRT
simpleRT_mean = squeeze(mean(simpleRT_data, 3, 'omitnan'));
simpleRT_se = squeeze(std(simpleRT_data, [], 3, 'omitnan'))/sqrt(size(simpleRT_data, 3));

% HEP gng
gng_mean = squeeze(mean(gng_data, 3, 'omitnan'));
gng_se = squeeze(std(gng_data, [], 3, 'omitnan'))/sqrt(size(gng_data, 3));

chan_dir = 'E:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study';
channels = load([chan_dir filesep 'chanlocs.mat']);
channels_num = [7, 16, 25, 34, 51];

for p = channels_num

    x_axis=51:2:400; % in mseconds
    figure;
    plot(x_axis, simpleRT_mean(p, :), 'color', [0 0.4470 0.7410],  'LineWidth',1.5) %plot mean gng (light blue)
    hold on
    jbfill(x_axis, [simpleRT_mean(p, :) + simpleRT_se(p, :)], [simpleRT_mean(p, :) - simpleRT_se(p, :)], [0 0.4470 0.7410],[0 0.4470 0.7410], 1, 0.2)
    hold on
    plot(x_axis, gng_mean(p, :), 'color', 'k',  'LineWidth',1.5) %plot mean gng (light blue)
    hold on
    jbfill(x_axis, [gng_mean(p, :) + gng_se(p, :)], [gng_mean(p, :) - gng_se(p, :)], 'k', 'k', 1, 0.2)
    hold on
    plot(x_axis, zeros(length(x_axis), 1), ':k')
    hold off
    box off
    ax = gca;
    if t == 1
        axis([-inf inf -1.5 1.5])
    else
        axis([-inf inf -2.5 1.5])
    end
    ax.LineWidth = 2.5;
    ax.FontSize = 28;
    ax.FontName = 'Arial';
    xlabel('Time (ms)', 'FontSize', 32, 'FontWeight','normal')
    ylabel('Amplitude (\muV)', 'FontSize', 32, 'FontWeight','normal')
    name_channels =  channels.chanlocs(p).labels;
    title(name_channels, 'FontSize', 32, 'FontWeight','normal')
%     legend(' before cue', 'b4', 'after cue', 'after')
end

%% plot young vs older - both tasks together - b4cue and after cue together
% tasks together
% y = channels x time x participants x measures
young_data = []; older_data =[];
young_data = squeeze(mean(y(:, :, gp == 1, :), 4));
older_data = squeeze(mean(y(:, :, gp == 2, :), 4));

% HEP young
young_mean = squeeze(mean(young_data, 3, 'omitnan'));
young_se = squeeze(std(young_data, [], 3, 'omitnan'))/sqrt(size(young_data, 3));

% HEP older
older_mean = squeeze(mean(older_data, 3, 'omitnan'));
older_se = squeeze(std(older_data, [], 3, 'omitnan'))/sqrt(size(older_data, 3));

chan_dir = 'E:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study';
channels = load([chan_dir filesep 'chanlocs.mat']);
channels_num = [7, 16, 25, 34, 51];

for p = channels_num

    x_axis=51:2:400; % in mseconds
    figure;
    plot(x_axis, young_mean(p, :), 'color', [0 0.4470 0.7410],  'LineWidth',1.5) %plot mean gng (light blue)
    hold on
    jbfill(x_axis, [young_mean(p, :) + young_se(p, :)], [young_mean(p, :) - young_se(p, :)], [0 0.4470 0.7410],[0 0.4470 0.7410], 1, 0.2)
    hold on
    plot(x_axis, older_mean(p, :), 'color', 'k',  'LineWidth',1.5) %plot mean gng (light blue)
    hold on
    jbfill(x_axis, [older_mean(p, :) + older_se(p, :)], [older_mean(p, :) - older_se(p, :)], 'k', 'k', 1, 0.2)
    hold on
    plot(x_axis, zeros(length(x_axis), 1), ':k')
    hold off
    box off
    ax = gca;
    if t == 1
        axis([-inf inf -1.5 1.5])
    else
        axis([-inf inf -2.5 1.5])
    end
    ax.LineWidth = 2.5;
    ax.FontSize = 28;
    ax.FontName = 'Arial';
    xlabel('Time (ms)', 'FontSize', 32, 'FontWeight','normal')
    ylabel('Amplitude (\muV)', 'FontSize', 32, 'FontWeight','normal')
    name_channels =  channels.chanlocs(p).labels;
    title(name_channels, 'FontSize', 32, 'FontWeight','normal')
%     legend(' before cue', 'b4', 'after cue', 'after')
end


%% plot data
% effect of before/after
figure
for chn=1:59
    plot(mean(squeeze(mean(y(chn, :, :, 1:2), 4)), 2))
    hold on
    plot(mean(squeeze(mean(y(chn, :, :, 3:4), 4)), 2))
    hold off
    legend('before', 'after')
    waitforbuttonpress
end
