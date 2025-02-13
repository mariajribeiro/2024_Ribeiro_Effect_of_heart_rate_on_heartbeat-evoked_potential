% calculate HEP for different heart rates - 25Jan2023 MariaRibeiro
clear; close all;
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
older = [11 12 14 20 21 22 23 32 37 38 41 43 47 48 49 52 55 57 58 63 64 65 67 69 7 70 71 75 8 83 86];
young = [13 15 16 25 26 28 31 33 34 36 4 42 44 45 46 50 51 53 54 56 59 6 62 66 68 72 74 76 78 80 82 84 85 9];

% all files
files_folder = 'E:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\eeg_derivatives';
folder_content = dir(files_folder);

% subjects included
subjects = {};
for f = 1:length(folder_content)
    if contains(folder_content(f).name, 'eeg_Rpeak_Tpeak')
        subjects = [subjects, folder_content(f).name(21:24)];
    end
end
subjects = unique(subjects);

tasks = {'passive', 'simpleRT', 'gonogo'};

HEP_IBI_tertiles = cell(2, 2); IBI_all = cell(2, 2);  IBI_tertiles = cell(2, 2); 
IBI_avg_std_min_max_range = cell(2, 2); %IBI_avg_std = cell(2, 2);
HR_avg_std_min_max_range = cell(2, 2);
yng = 0; old = 0;
for s = 1:length(subjects)
    if strcmp(subjects{s}(4), '_')
        sbj = 3;
    else
        sbj = 3:4;
    end
    
    subj_number = str2num(subjects{s}(sbj));
    if ismember(subj_number, young)
        group = 1; yng = yng + 1;
    else
        group = 2; old = old + 1;
    end
    
    HEP_IBI_tertiles{group, 1} = [ HEP_IBI_tertiles{group, 1}; subj_number];
    IBI_all{group, 1} = [ IBI_all{group, 1}; subj_number]; 
%     IBI_avg_std{group, 1} = [ IBI_avg_std{group, 1}; subj_number]; 
    IBI_tertiles{group, 1} = [IBI_tertiles{group, 1}; subj_number]; 
    IBI_avg_std_min_max_range{group, 1} = [IBI_avg_std_min_max_range{group, 1}; subj_number]; 
    HR_avg_std_min_max_range{group, 1} = [HR_avg_std_min_max_range{group, 1}; subj_number]; 
       
    for t = 1:length(tasks)
        % clear eeglab
        STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];  

        % load data
        filename = {};
        for f = 1:length(folder_content)
            if contains(folder_content(f).name, subjects(s)) && contains(folder_content(f).name, tasks(t))
                filename = folder_content(f).name;
                EEG = pop_loadset('filename', filename, 'filepath', files_folder );
%                 EEG = pop_select( EEG, 'nochannel',{'EKG'});
                [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
            end
        end

        % append data from both runs
        if length(ALLEEG) == 2
            EEG = pop_mergeset( ALLEEG, [1  2], 0);
            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'gui','off'); 
        end
        
        % delete R and T peaks associated with abnormal cardiac cycles or
        % cycles where the cycle was not segmented properly
        outliers_evnt = find_cardiac_cycle_outliers(EEG);
        % delete events associated with outlier cycles
        EEG = pop_selectevent( EEG, 'omitevent', outliers_evnt ,'deleteevents','on');
        
        % make epoch with all R-peaks
        EEG = pop_epoch( EEG, {  'Rpeak'  }, [-0.15         1.5], 'newname', 'Rpeak epochs all', 'epochinfo', 'yes');
        EEG = pop_rmbase( EEG, [-150 -50] ,[]);
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
        
        % calculate heart rate for each epoch
        IBI = [];
        for ep = 1:length(EEG.epoch)
            Rpeaks = find(strcmp(EEG.epoch(ep).eventtype, 'Rpeak'));
            if length(Rpeaks) > 1 && EEG.epoch(ep).eventlatency{Rpeaks(2)} > 0
                IBI(ep) = EEG.epoch(ep).eventlatency{Rpeaks(2)};
            else
                IBI(ep) = NaN;
            end
        end
                
        % find outliers
        zscore_IBI = (IBI - mean(IBI, 'omitnan'))/std(IBI, 'omitnan');
        IBI_outliers = find(abs(zscore_IBI) > 4);
        IBI(IBI_outliers) = NaN;
        
        % delete outliers
        HEP_data = EEG.data;
        HEP_data(:, :, isnan(IBI)) = [];
        IBI(isnan(IBI)) = [];

        % divide HEP according to heart rate (interbeat interval) tertiles
        [IBI_sorted,index] = sortrows(IBI');
        
        IBI1 = mean(IBI(index(1:floor(length(index)/3))));
        IBI2 = mean(IBI(index(floor(length(index)/3)+1:2*floor(length(index)/3))));
        IBI3 = mean(IBI(index(2*floor(length(index)/3)+1:end)));
        
        HEP1 = mean(HEP_data(:, :, index(1:floor(length(index)/3))), 3);
        HEP2 = mean(HEP_data(:, :, index(floor(length(index)/3)+1:2*floor(length(index)/3))), 3);
        HEP3 = mean(HEP_data(:, :, index(2*floor(length(index)/3)+1:end)), 3);
        
        HR = IBI.^-1*60000;
        
        if group == 1
            IBI_tertiles{group, 2}(yng, t, 1) = IBI1;
            IBI_tertiles{group, 2}(yng, t, 2) = IBI2;
            IBI_tertiles{group, 2}(yng, t, 3) = IBI3;
            
            HEP_IBI_tertiles{group, 2}(yng, t, 1, :, :) = HEP1;
            HEP_IBI_tertiles{group, 2}(yng, t, 2, :, :) = HEP2;
            HEP_IBI_tertiles{group, 2}(yng, t, 3, :, :) = HEP3;
            
            IBI_all{group, 2}{yng, t} = IBI;
            IBI_avg_std_min_max_range{group, 2}(yng, t, :) = [mean(IBI), std(IBI), min(IBI), max(IBI), max(IBI)-min(IBI)];
            HR_avg_std_min_max_range{group, 2}(yng, t, :) = [mean(HR), std(HR), min(HR), max(HR), max(HR)-min(HR)];
            
        else
            IBI_tertiles{group, 2}(old, t, 1) = IBI1;
            IBI_tertiles{group, 2}(old, t, 2) = IBI2;
            IBI_tertiles{group, 2}(old, t, 3) = IBI3;
            
            HEP_IBI_tertiles{group, 2}(old, t, 1, :, :) = HEP1;
            HEP_IBI_tertiles{group, 2}(old, t, 2, :, :) = HEP2;
            HEP_IBI_tertiles{group, 2}(old, t, 3, :, :) = HEP3;
            
            IBI_all{group, 2}{old, t}= IBI;
            IBI_avg_std_min_max_range{group, 2}(old, t, :) = [mean(IBI), std(IBI), min(IBI), max(IBI), max(IBI)-min(IBI)];
            HR_avg_std_min_max_range{group, 2}(old, t, :) = [mean(HR), std(HR), min(HR), max(HR), max(HR)-min(HR)];
        end
    end       
end

save HEP_IBI_tertiles HEP_IBI_tertiles
save IBI_tertiles IBI_tertiles
save IBI_all IBI_all
save IBI_avg_std_min_max_range IBI_avg_std_min_max_range
save HR_avg_std_min_max_range HR_avg_std_min_max_range
%% plot figures
% plot average heart rate both tasks 2 groups - simple RT vs gng
load IBI_all; load IBI_avg_std_min_max_range
% plot_2conditions(grp1_cond1, grp1_cond2, grp2_cond1, grp2_cond2, y_label_text, x_label, title_text)
plot_2conditions(IBI_avg_std_min_max_range{1, 2}(:,2, 1), IBI_avg_std_min_max_range{1, 2}(:,3, 1),...
    IBI_avg_std_min_max_range{2, 2}(:,2, 1), IBI_avg_std_min_max_range{2, 2}(:,3, 1), 'IBI (ms)', {'Simple RT' 'Go/no-go' 'Simple RT' 'Go/no-go'},...
    '')

plot_2conditions(IBI_avg_std_min_max_range{1, 2}(:,2, 2), IBI_avg_std_min_max_range{1, 2}(:,3, 2),...
    IBI_avg_std_min_max_range{2, 2}(:,2, 2), IBI_avg_std_min_max_range{2, 2}(:,3, 2), 'IBI std (ms)', {'Simple RT' 'Go/no-go' 'Simple RT' 'Go/no-go'},...
    '')

plot_2conditions(IBI_avg_std_min_max_range{1, 2}(:,2, 5), IBI_avg_std_min_max_range{1, 2}(:,3, 5),...
    IBI_avg_std_min_max_range{2, 2}(:,2, 5), IBI_avg_std_min_max_range{2, 2}(:,3, 5), 'IBI range (ms)', {'Simple RT' 'Go/no-go' 'Simple RT' 'Go/no-go'},...
    '')

%% plot as heart rate not IBI

% plot_2conditions(grp1_cond1, grp1_cond2, grp2_cond1, grp2_cond2, y_label_text, x_label, title_text)
plot_2conditions(HR_avg_std_min_max_range{1, 2}(:,2, 1), HR_avg_std_min_max_range{1, 2}(:,3, 1),...
    HR_avg_std_min_max_range{2, 2}(:,2, 1), HR_avg_std_min_max_range{2, 2}(:,3, 1), 'HR (bpm)', {'Simple RT' 'Go/no-go' 'Simple RT' 'Go/no-go'},...
    '')

plot_2conditions(HR_avg_std_min_max_range{1, 2}(:,2, 2), HR_avg_std_min_max_range{1, 2}(:,3, 2),...
    HR_avg_std_min_max_range{2, 2}(:,2, 2), HR_avg_std_min_max_range{2, 2}(:,3, 2), 'HR std (bpm)', {'Simple RT' 'Go/no-go' 'Simple RT' 'Go/no-go'},...
    '')

plot_2conditions(HR_avg_std_min_max_range{1, 2}(:,2, 5), HR_avg_std_min_max_range{1, 2}(:,3, 5),...
    HR_avg_std_min_max_range{2, 2}(:,2, 5), HR_avg_std_min_max_range{2, 2}(:,3, 5), 'HR range (bpm)', {'Simple RT' 'Go/no-go' 'Simple RT' 'Go/no-go'},...
    '')

%% save HR data in excel file to perform stats in SPSS
load HR_avg_std_min_max_range
% create table
% participant, group, IBI_avg_passive, IBI_avg_simpleRT, IBI_avg_gng, IBI_std_passive, IBI_std_simpleRT, IBI_std_gng
% IBI_avg_std{group, 2}(old, t, :) = [mean(IBI), std(IBI)];
HR_table = table([HR_avg_std_min_max_range{1, 1};HR_avg_std_min_max_range{2, 1}],...
    [ones(length(HR_avg_std_min_max_range{1, 1}), 1); ones(length(HR_avg_std_min_max_range{2, 1}), 1)*2], ...
    [HR_avg_std_min_max_range{1, 2}(:, 1, 1); HR_avg_std_min_max_range{2, 2}(:, 1, 1)], ...
    [HR_avg_std_min_max_range{1, 2}(:, 2, 1); HR_avg_std_min_max_range{2, 2}(:, 2, 1)], ...
    [HR_avg_std_min_max_range{1, 2}(:, 3, 1); HR_avg_std_min_max_range{2, 2}(:, 3, 1)], ...
    [HR_avg_std_min_max_range{1, 2}(:, 1, 2); HR_avg_std_min_max_range{2, 2}(:, 1, 2)], ...
    [HR_avg_std_min_max_range{1, 2}(:, 2, 2); HR_avg_std_min_max_range{2, 2}(:, 2, 2)], ...
    [HR_avg_std_min_max_range{1, 2}(:, 3, 2); HR_avg_std_min_max_range{2, 2}(:, 3, 2)], ...
    [HR_avg_std_min_max_range{1, 2}(:, 1, 5); HR_avg_std_min_max_range{2, 2}(:, 1, 5)], ...
    [HR_avg_std_min_max_range{1, 2}(:, 2, 5); HR_avg_std_min_max_range{2, 2}(:, 2, 5)], ...
    [HR_avg_std_min_max_range{1, 2}(:, 3, 5); HR_avg_std_min_max_range{2, 2}(:, 3, 5)], ...
    'VariableNames',{'participant', 'group', 'HR_avg_passive', 'HR_avg_simpleRT', 'HR_avg_gng',...
    'HR_std_passive', 'HR_std_simpleRT', 'HR_std_gng', ...
    'HR_range_passive', 'HR_range_simpleRT', 'HR_range_gng'});

filename = 'HR_avg_std_range.xlsx';
writetable(HR_table,filename)

%% plot for each participant IBI across cardiac cycles for simple RT run 1
load IBI_all; load IBI_avg_std
% IBI_avg_std{group, 2}(old, t, :)
for grp = 1:2
    for p = 1:length(IBI_avg_std{1, 1})
        clear x
    %     IBI_all{group, 2}{old, t}= IBI;
        IBI = IBI_all{grp, 2}{p, 2};
        x(1) = IBI(1);
        for i =1:length(IBI)-1
            x(i+1) = x(i)+IBI(i+1);
        end

        figure;
        plot(x/1000, IBI, '-o', 'color', [.5 .5 .5], ...
            'MarkerEdgeColor','k','MarkerFaceColor','k', 'MarkerSize',2, 'LineWidth', .5);
        axis([0 1000 550 1000]);
        ax = gca; box off
        ax.LineWidth = 2.5;
        ax.FontSize = 24;
        ax.FontName = 'Arial';
        ax.Color = 'none';
%         ax.XAxis.FontSize = 20;
        xlabel('Time (s)', 'FontSize', 30, 'FontWeight','normal')
        ylabel('IBI (ms)', 'FontSize', 30, 'FontWeight','normal')
        title(['AB', num2str(IBI_avg_std{grp, 1}(p)), 'group', num2str(grp)], 'FontSize', 32, 'FontWeight','normal')
        waitforbuttonpress
    end
end


%% save data in excel file to perform stats in SPSS
load IBI_avg_std
% create table
% participant, group, IBI_avg_passive, IBI_avg_simpleRT, IBI_avg_gng, IBI_std_passive, IBI_std_simpleRT, IBI_std_gng
% IBI_avg_std{group, 2}(old, t, :) = [mean(IBI), std(IBI)];
IBI_table = table([IBI_avg_std{1, 1};IBI_avg_std{2, 1}],...
    [ones(length(IBI_avg_std{1, 1}), 1); ones(length(IBI_avg_std{2, 1}), 1)*2], ...
    [IBI_avg_std{1, 2}(:, 1, 1); IBI_avg_std{2, 2}(:, 1, 1)], ...
    [IBI_avg_std{1, 2}(:, 2, 1); IBI_avg_std{2, 2}(:, 2, 1)], ...
    [IBI_avg_std{1, 2}(:, 3, 1); IBI_avg_std{2, 2}(:, 3, 1)], ...
    [IBI_avg_std{1, 2}(:, 1, 2); IBI_avg_std{2, 2}(:, 1, 2)], ...
    [IBI_avg_std{1, 2}(:, 2, 2); IBI_avg_std{2, 2}(:, 2, 2)], ...
    [IBI_avg_std{1, 2}(:, 3, 2); IBI_avg_std{2, 2}(:, 3, 2)], ...
    'VariableNames',{'participant', 'group', 'IBI_avg_passive', 'IBI_avg_simpleRT', 'IBI_avg_gng',...
    'IBI_std_passive', 'IBI_std_simpleRT', 'IBI_std_gng'});

filename = 'IBI_avg_std.xlsx';
writetable(IBI_table,filename)


load IBI_avg_std_min_max_range
% create table
% participant, group, IBI_avg_passive, IBI_avg_simpleRT, IBI_avg_gng,
% IBI_std_passive, IBI_std_simpleRT, IBI_std_gng,
% IBI_min_passive, IBI_min_simpleRT, IBI_min_gng
% IBI_max_passive, IBI_max_simpleRT, IBI_max_gng
% IBI_range_passive, IBI_range_simpleRT, IBI_range_gng

IBI_table = table([IBI_avg_std_min_max_range{1, 1};IBI_avg_std_min_max_range{2, 1}],...
    [ones(length(IBI_avg_std_min_max_range{1, 1}), 1); ones(length(IBI_avg_std_min_max_range{2, 1}), 1)*2], ...
    [IBI_avg_std_min_max_range{1, 2}(:, 1, 1); IBI_avg_std_min_max_range{2, 2}(:, 1, 1)], ...
    [IBI_avg_std_min_max_range{1, 2}(:, 2, 1); IBI_avg_std_min_max_range{2, 2}(:, 2, 1)], ...
    [IBI_avg_std_min_max_range{1, 2}(:, 3, 1); IBI_avg_std_min_max_range{2, 2}(:, 3, 1)], ...
    [IBI_avg_std_min_max_range{1, 2}(:, 1, 2); IBI_avg_std_min_max_range{2, 2}(:, 1, 2)], ...
    [IBI_avg_std_min_max_range{1, 2}(:, 2, 2); IBI_avg_std_min_max_range{2, 2}(:, 2, 2)], ...
    [IBI_avg_std_min_max_range{1, 2}(:, 3, 2); IBI_avg_std_min_max_range{2, 2}(:, 3, 2)], ...
        [IBI_avg_std_min_max_range{1, 2}(:, 1, 3); IBI_avg_std_min_max_range{2, 2}(:, 1, 3)], ...
    [IBI_avg_std_min_max_range{1, 2}(:, 2, 3); IBI_avg_std_min_max_range{2, 2}(:, 2, 3)], ...
    [IBI_avg_std_min_max_range{1, 2}(:, 3, 3); IBI_avg_std_min_max_range{2, 2}(:, 3, 3)], ...
        [IBI_avg_std_min_max_range{1, 2}(:, 1, 4); IBI_avg_std_min_max_range{2, 2}(:, 1, 4)], ...
    [IBI_avg_std_min_max_range{1, 2}(:, 2, 4); IBI_avg_std_min_max_range{2, 2}(:, 2, 4)], ...
    [IBI_avg_std_min_max_range{1, 2}(:, 3, 4); IBI_avg_std_min_max_range{2, 2}(:, 3, 4)], ...
        [IBI_avg_std_min_max_range{1, 2}(:, 1, 5); IBI_avg_std_min_max_range{2, 2}(:, 1, 5)], ...
    [IBI_avg_std_min_max_range{1, 2}(:, 2, 5); IBI_avg_std_min_max_range{2, 2}(:, 2, 5)], ...
    [IBI_avg_std_min_max_range{1, 2}(:, 3, 5); IBI_avg_std_min_max_range{2, 2}(:, 3, 5)], ...
    'VariableNames',{'participant', 'group', 'IBI_avg_passive', 'IBI_avg_simpleRT', 'IBI_avg_gng',...
    'IBI_std_passive', 'IBI_std_simpleRT', 'IBI_std_gng', ...
    'IBI_min_passive', 'IBI_min_simpleRT', 'IBI_min_gng',...
    'IBI_max_passive', 'IBI_max_simpleRT', 'IBI_max_gng',...
    'IBI_range_passive', 'IBI_range_simpleRT', 'IBI_range_gng'});

filename = 'IBI_avg_std_min_max_range.xlsx';
writetable(IBI_table,filename)

%%
load IBI_tertiles
colormap cool;
cmap = colormap;
colors_ibi = [cmap(1,:); cmap(128,:); cmap(256,:)].*255;

tasks = {'passive', 'simple RT', 'go/no-go'};
for t = 1:length(tasks)
   figure;
   mean_ibi = squeeze(mean(IBI_tertiles{1, 2}(:, t, :), 1));
   se_ibi = squeeze(std(IBI_tertiles{1, 2}(:, t, :), [], 1))./sqrt(size(IBI_tertiles{1, 2}, 1));
   errorbar([1:3], mean_ibi, se_ibi, '-o', 'color', 'k', 'LineWidth', 2)
   hold on 
   mean_ibi = squeeze(mean(IBI_tertiles{2, 2}(:, t, :), 1));
   se_ibi = squeeze(std(IBI_tertiles{2, 2}(:, t, :), [], 1))./sqrt(size(IBI_tertiles{1, 2}, 1));
   errorbar([1:3], mean_ibi, se_ibi, '-d', 'color', 'r', 'LineWidth', 2)
    hold off
    box off
    ax = gca;
    ax.LineWidth = 2.5;
    ax.FontSize = 28;
    ax.FontName = 'Arial';
   title([tasks{t}] , 'FontSize', 32, 'FontWeight','normal')
   axis([0 4 -inf inf])
   ylabel('Interbeat interval (ms)', 'FontSize', 24, 'FontWeight','normal')
   set(gca, 'XTick',1:3, 'XTickLabel', 1:3)
   set(gca, 'XLim',[0 4])
end



%%
channels = load('chanlocs_with_ekg.mat');
tasks = {'passive', 'simple RT', 'go/nogo'};
load HEP_IBI_tertiles;

x_axis=-0.148:0.002:1.5; % in seconds
% time window where ECG shows a significant effect of IBI
mask_task{2} = find(x_axis <= .532);
mask_task{3} = find(x_axis <= .476);

% load time windows where HEP correlates with IBI
% masks were calculated from 100-700 ms after R peak
mask_dir = 'E:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\LIMO_stats\HEPvsIBI\level2_onesample_IBIvsHEP';
mask_simpleRT = load([mask_dir, '\simpleRT\mask']);
mask_gng = load([mask_dir, '\gng\mask']);
mask_ecg_dir = 'E:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\LIMO_stats\ECGvsIBI\level2_onesample_IBIvsECG';
mask_ecg_simpleRT = load([mask_ecg_dir, '\simpleRT\mask']);
mask_ecg_gng = load([mask_ecg_dir, '\gng\mask']);

for t = 1:length(tasks)
% 
% figure;
% for p = 1:size(HEP_b4cue{2, 2}, 1)
%     plot(squeeze(HEP_b4cue{2, 2}(p, t, 7, :)))
%     waitforbuttonpress
% end
% close all

for ibi = 1:3
    young_ibi_mean(ibi, :, :) = squeeze(mean(HEP_IBI_tertiles{1, 2}(:, t, ibi, :, :), 1));
    older_ibi_mean(ibi, :, :) = squeeze(mean(HEP_IBI_tertiles{2, 2}(:, t, ibi, :, :), 1));

    young_ibi_se(ibi, :, :) = squeeze(std(HEP_IBI_tertiles{1, 2}(:, t, ibi, :, :), [], 1))/sqrt(size(HEP_IBI_tertiles{1, 2}, 1));
    older_ibi_se(ibi, :, :)= squeeze(std(HEP_IBI_tertiles{2, 2}(:, t, ibi, :, :), [], 1))/sqrt(size(HEP_IBI_tertiles{1, 2}, 1));
end
% figure;
% topoplot(young_b4cue_mean(:, 225), channels.chanlocs); 
% Axis.FontSize = 16; caxis([-2 2]); 
% colorbar;
% % colorbar('Ticks',[0, 2, 4, 6], 'FontSize', 30, 'FontWeight','normal');
% colormap(crameri('hawaii'));
% title(['young b4cue ', tasks{t}] , 'FontSize', 32, 'FontWeight','normal')
% 
% figure;
% topoplot(older_b4cue_mean(:, 225), channels.chanlocs); 
% Axis.FontSize = 16; caxis([-2 2]); 
% colorbar;
% % colorbar('Ticks',[0, 2, 4, 6], 'FontSize', 30, 'FontWeight','normal');
% colormap(crameri('hawaii'));
% title(['older b4cue ', tasks{t}] , 'FontSize', 32, 'FontWeight','normal')



    % channels = F3,Fz,F4,C3,Cz,C4,P3,Pz,P4
    % channels_num = [7, 16, 25, 34];
    channels_num = 60;%34;%[7, 25, 43, 57, 60];
    colormap cool;
    cmap = colormap;
    colors_ibi = [cmap(1,:); cmap(size(cmap, 1)/2,:); cmap(end,:)];
    for p = channels_num
        name_channels = channels.chanlocs_with_ekg(p).labels;

        figure;
        for ibi = 1:3
            plot(x_axis, squeeze(young_ibi_mean(ibi, p, :)), 'color', colors_ibi(ibi, :),  'LineWidth',2.5)
            hold on
            jbfill(x_axis, squeeze([young_ibi_mean(ibi, p, :) + young_ibi_se(ibi, p, :)])', squeeze([young_ibi_mean(ibi, p, :) - young_ibi_se(ibi, p, :)])', colors_ibi(ibi, :), colors_ibi(ibi, :), 1, .2)
            hold on
        end
        plot(x_axis, zeros(length(x_axis), 1), '--k')
        if strcmp(name_channels, 'EKG')
            mask = [];
            if t == 2
                mask = find(mask_ecg_simpleRT.mask == 1);
            elseif t == 3
                mask = find(mask_ecg_simpleRT.mask == 1);
            end
                plot(x_axis(mask+125), ones(length(mask), 1)  * -50, 's', 'MarkerSize', 4.5, 'color', [.5 .5 .5], ...
            'MarkerFaceColor', [.5 .5 .5])
        elseif t == 2 % simple RT
            mask = find(mask_simpleRT.mask(p, :) == 1);
             plot(x_axis(mask+125), ones(length(mask), 1) * -1, 's', 'MarkerSize', 4.5, 'color', [.5 .5 .5], ...
            'MarkerFaceColor', [.5 .5 .5])
        elseif t == 3 % gng
            mask = find(mask_gng.mask(p, :) == 1);
             plot(x_axis(mask+125), ones(length(mask), 1) * -1, 's', 'MarkerSize', 4.5, 'color', [.5 .5 .5], ...
            'MarkerFaceColor', [.5 .5 .5])
        end
        hold off
        box off
        ax = gca;
        if strcmp(name_channels, 'EKG')
            axis([-inf .7 -100 250])%axis([-inf .7 -200 200])
        else
            axis([-inf .7 -1.50 1.50]);%axis([-inf .7 -1.50 1.50])
        end
%         if t == 1
%             axis([-inf .7 -1.5 1])
%         else
%             axis([-inf .7 -1.5 1])
%         end
        ax.LineWidth = 2.5;
        ax.FontSize = 28;
        ax.FontName = 'Arial';
        xlabel('Time (s)', 'FontSize', 32, 'FontWeight','normal')
        ylabel('Amplitude (\muV)', 'FontSize', 32, 'FontWeight','normal')
        name_channels =  channels.chanlocs_with_ekg(p).labels;
        title(['Young - ', tasks{t}, ' - ', name_channels] , 'FontSize', 28, 'FontWeight','normal')
%         basefilename2 = sprintf('Tepochplot_gng-sRT_young_cut_%d.fig', p);
%         savefig(basefilename2);


        figure;
        for ibi = 1:3
            plot(x_axis, squeeze(older_ibi_mean(ibi, p, :)), 'color', colors_ibi(ibi, :),  'LineWidth',2.5)
            hold on
            jbfill(x_axis, squeeze([older_ibi_mean(ibi, p, :) + older_ibi_se(ibi, p, :)])', squeeze([older_ibi_mean(ibi, p, :) - older_ibi_se(ibi, p, :)])', colors_ibi(ibi, :), colors_ibi(ibi, :), 1, .2)
            hold on
        end
        plot(x_axis, zeros(length(x_axis), 1), '--k')
        if strcmp(name_channels, 'EKG')
            mask = [];
            if t == 2
                mask = find(mask_ecg_simpleRT.mask == 1);
            elseif t == 3
                mask = find(mask_ecg_simpleRT.mask == 1);
            end
                plot(x_axis(mask+125), ones(length(mask), 1)  * -50, 's', 'MarkerSize', 4.5, 'color', [.5 .5 .5], ...
            'MarkerFaceColor', [.5 .5 .5])
        elseif t == 2 % simple RT
            mask = find(mask_simpleRT.mask(p, :) == 1);
            plot(x_axis(mask+125), ones(length(mask), 1) * -1, 's', 'MarkerSize', 4.5, 'color', [.5 .5 .5], ...
            'MarkerFaceColor', [.5 .5 .5])
        elseif t == 3 % gng
            mask = find(mask_gng.mask(p, :) == 1);
             plot(x_axis(mask+125), ones(length(mask), 1) * -1, 's', 'MarkerSize', 4.5, 'color', [.5 .5 .5], ...
            'MarkerFaceColor', [.5 .5 .5])
        end
        hold off
        box off
        ax = gca;
        if p == 60 
            axis([-inf .7  -100 250])%axis([-inf .7 -200 200])
        else
            axis([-inf .7  -1.50 1.50])
        end
%         if t == 1
%             axis([-inf inf -1.5 1.5])
%         else
%             axis([-inf inf -2.5 1.5])
%         end
        ax.LineWidth = 2.5;
        ax.FontSize = 28;
        ax.FontName = 'Arial';
        xlabel('Time (s)', 'FontSize', 32, 'FontWeight','normal')
        ylabel('Amplitude (\muV)', 'FontSize', 32, 'FontWeight','normal')
        name_channels =  channels.chanlocs_with_ekg(p).labels;
        title(['Older - ', tasks{t}, ' - ', name_channels] , 'FontSize', 28, 'FontWeight','normal')

    end

end

%% plot effect of IBI on ECG - average of both groups, average of both tasks

load HEP_IBI_tertiles;



% load time windows where HEP correlates with IBI
% masks were calculated from 100-700 ms after R peak
mask_ecg_dir = 'E:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\LIMO_stats\ECGvsIBI\level2_onesample_IBIvsECG';
mask_ecg_simpleRT = load([mask_ecg_dir, '\simpleRT\mask']);
mask_ecg_gng = load([mask_ecg_dir, '\gng\mask']);
%     HEP_IBI_tertiles = participants x tasks x ibi x channels x time
%avg data across tasks
young_data = squeeze(mean(HEP_IBI_tertiles{1, 2}(:, 2:3, :, 60, :), 2));
older_data = squeeze(mean(HEP_IBI_tertiles{2, 2}(:, 2:3, :, 60, :), 2));

data_all = cat(1, young_data,  older_data);

data_mean = squeeze(mean(data_all, 1));
data_se = squeeze(std(data_all, [], 1))/sqrt(size(data_all, 1));

colormap cool;
cmap = colormap;
colors_ibi = [cmap(1,:); cmap(size(cmap, 1)/2,:); cmap(end,:)];
x_axis=-0.148:0.002:1.5; % in seconds
figure;
for ibi = 1:3
    plot(x_axis, squeeze(data_mean(ibi, :)), 'color', colors_ibi(ibi, :),  'LineWidth',2.5)
    hold on
    jbfill(x_axis, squeeze(data_mean(ibi, :) + data_se(ibi, :)), squeeze(data_mean(ibi, :) - data_se(ibi, :)), colors_ibi(ibi, :), colors_ibi(ibi, :), 1, .2)
    hold on
end
plot(x_axis, zeros(length(x_axis), 1), '--k')
mask = [];
mask = find(mask_ecg_simpleRT.mask == 1);
plot(x_axis(mask+125), ones(length(mask), 1)  * -50, 's', 'MarkerSize', 4.5, 'color', [.5 .5 .5], ...
'MarkerFaceColor', [.5 .5 .5])
hold on
mask = find(mask_ecg_simpleRT.mask == 1);
plot(x_axis(mask+125), ones(length(mask), 1)  * -50, 's', 'MarkerSize', 4.5, 'color', [.5 .5 .5], ...
'MarkerFaceColor', [.5 .5 .5])

hold off
box off
ax = gca;
    axis([-inf .7 -250 250]) % axis([-inf .7 -200 200])
ax.LineWidth = 2.5;
ax.FontSize = 28;
ax.FontName = 'Arial';
xlabel('Time (s)', 'FontSize', 32, 'FontWeight','normal')
ylabel('Amplitude (\muV)', 'FontSize', 32, 'FontWeight','normal')
title('ECG', 'FontSize', 28, 'FontWeight','normal')


%%
function plot_2conditions(grp1_cond1, grp1_cond2, grp2_cond1, grp2_cond2, y_label_text, x_label, title_text)
    figure;
    index_nan1 = find(isnan(grp1_cond1));
    index_nan2 = find(isnan(grp1_cond2));

    grp1_cond1(unique([index_nan1, index_nan2])) = [];
    grp1_cond2(unique([index_nan1, index_nan2])) = [];
    
    % individual points
    for y=1:length(grp1_cond1)
    plot([1+rand*0.2-0.1, 2+rand*0.2-0.1], [grp1_cond1(y), grp1_cond2(y)],'-o', 'color', [.5 .5 .5], ...
        'MarkerEdgeColor',[.5 .5 .5],'MarkerSize',8, 'LineWidth', 1.5);
    hold on;
    end

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
    

    
    % plot data for group 2
    index_nan1 = find(isnan(grp2_cond1));
    index_nan2 = find(isnan(grp2_cond2));

    grp2_cond1(unique([index_nan1, index_nan2])) = [];
    grp2_cond2(unique([index_nan1, index_nan2])) = [];
    
    % individual points
    for y=1:length(grp2_cond1)
        plot([4+rand*0.2-0.1, 5+rand*0.2-0.1], [grp2_cond1(y), grp2_cond2(y)],'-o', 'color', [.5 .5 .5], ...
            'MarkerEdgeColor', [.5 .5 .5],'MarkerSize',8, 'LineWidth', 1.5);
        hold on;
    end

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