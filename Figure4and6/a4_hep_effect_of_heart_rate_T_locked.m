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

% subjects excluded from T-locked HEP analyses because T waves were not
% segmented properly = AB21, AB50 and AB55
excluded_sbjts = {'AB21', 'AB50', 'AB55'};


tasks = {'passive', 'simpleRT', 'gonogo'};

HEP_T_locked_IBI_tertiles = cell(2, 2);
yng = 0; old = 0;
for s = 1:length(subjects)
    if ismember(subjects{s}, excluded_sbjts)
        continue
    end
    if strcmp(subjects{s}(4), '_')
        sbj = 3;
    else
        sbj = 3:4;
    end
    
    subj_number = str2double(subjects{s}(sbj));
    if ismember(subj_number, young)
        group = 1; yng = yng + 1;
    else
        group = 2; old = old + 1;
    end
    
    HEP_T_locked_IBI_tertiles{group, 1} = [ HEP_T_locked_IBI_tertiles{group, 1}; subj_number];
       
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
        
        % make epoch with all R-peaks
        EEG = pop_epoch( EEG, {  'Rpeak'  }, [-0.15         1.5], 'newname', 'Rpeak epochs all', 'epochinfo', 'yes');
        EEG = pop_rmbase( EEG, [-150 -50] ,[]);
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
        
         % calculate heart rate for each epoch
        IBI_tmp = []; RT_interval= [];
        for ep = 1:length(EEG.epoch)
            Rpeaks = find(strcmp(EEG.epoch(ep).eventtype, 'Rpeak'));
            Tpeaks = find(strcmp(EEG.epoch(ep).eventtype, 'Tpeak'));
            if length(Rpeaks) > 1 && length(Tpeaks) >= 1 && Tpeaks(1) < Rpeaks(2)
                IBI_tmp(ep) = EEG.epoch(ep).eventlatency{Rpeaks(2)};
                RT_interval(ep) = EEG.epoch(ep).eventlatency{Tpeaks(1)};
            else
                IBI_tmp(ep) = NaN;
                RT_interval(ep) = NaN;
            end
        end
                
        % find outliers - ectopic beats or cycles where ECG was not
        % properly segmented
        zscore_IBI = (IBI_tmp - mean(IBI_tmp, 'omitnan'))/std(IBI_tmp, 'omitnan');
        IBI_tmp(abs(zscore_IBI) > 4) = NaN;
        zscore_RT_interval = (RT_interval - mean(RT_interval, 'omitnan'))/std(RT_interval, 'omitnan');
        RT_interval(abs(zscore_RT_interval) > 4) = NaN;
        
        
        % Exclude outliers that might refer to ectopic beats or
        % trials where the Rpeak/Tpeak was wrongly detected
        trials2keep = intersect(find(~isnan(IBI_tmp)), find(~isnan(RT_interval)));
        
        EEG = pop_select( EEG, 'trial', trials2keep);
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 4,'gui','off'); 
        
        % exclude outliers from IBI variable
        IBI = IBI_tmp(trials2keep)';
        
        % in order to ensure that each epoch is associated with one
        % T-locked HEP and one IBI we need to remove teh second T event
        % from each R-locked epoch
        Tevent = []; 
        for evnt = 1:length(EEG.event)
            if strcmp(EEG.event(evnt).type,  'Tpeak')
                Tevent = [Tevent; evnt EEG.event(evnt).epoch];
            end
        end
        evnt2del =[];
        for evnt = 2:size(Tevent, 1)
             if Tevent(evnt, 2) == Tevent(evnt-1, 2)
                    evnt2del = [evnt2del; Tevent(evnt, 1)];
             end
        end
        EEG.event(evnt2del) = [];
        
        t_peaks = 0;
        for evn = 1:length(EEG.event)
            if strcmp(EEG.event(evn).type, 'Tpeak')
                t_peaks = t_peaks + 1;
            end
        end
        if t_peaks ~= 0 % in case t peaks were detected
            EEG = pop_epoch( EEG, {  'Tpeak'  }, [-0.2         0.4], 'newname', 'Tpeak epochs b4cue', 'epochinfo', 'yes');
            [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
            % HEP T-locked data
            HEP_T_locked_data = EEG.data(1:60, :, :);
            number_trials(t) = size(EEG.data, 3);
        else % in case t peaks were not detected
            % HEP T-locked data
            HEP_T_locked_data = [];
            number_trials(t) = 0;
        end

        % divide HEP according to heart rate (interbeat interval) tertiles
        [IBI_sorted,index] = sortrows(IBI);
        
        HEP1 = mean(HEP_T_locked_data(:, :, index(1:floor(length(index)/3))), 3);
        HEP2 = mean(HEP_T_locked_data(:, :, index(floor(length(index)/3)+1:2*floor(length(index)/3))), 3);
        HEP3 = mean(HEP_T_locked_data(:, :, index(2*floor(length(index)/3)+1:end)), 3);
        

        
        if group == 1
            HEP_T_locked_IBI_tertiles{group, 2}(yng, t, 1, :, :) = HEP1;
            HEP_T_locked_IBI_tertiles{group, 2}(yng, t, 2, :, :) = HEP2;
            HEP_T_locked_IBI_tertiles{group, 2}(yng, t, 3, :, :) = HEP3;
        else
            HEP_T_locked_IBI_tertiles{group, 2}(old, t, 1, :, :) = HEP1;
            HEP_T_locked_IBI_tertiles{group, 2}(old, t, 2, :, :) = HEP2;
            HEP_T_locked_IBI_tertiles{group, 2}(old, t, 3, :, :) = HEP3;
        end   
    end
end

save HEP_T_locked_IBI_tertiles HEP_T_locked_IBI_tertiles

%%
channels = load('chanlocs_with_ekg.mat');
tasks = {'passive', 'simple RT', 'go/no-go'};
load HEP_T_locked_IBI_tertiles;

x_axis=-199:2:400; % in miliseconds

% load time windows where HEP correlates with IBI
mask_dir = 'E:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\LIMO_stats\HEP_TLockedvsIBI\level2_onesample_IBIvsHEP';
mask_simpleRT = load([mask_dir, '\simpleRT\mask']);
mask_gng = load([mask_dir, '\gng\mask']);

for t = 1:length(tasks)

for ibi = 1:3
    young_ibi_mean(ibi, :, :) = squeeze(mean(HEP_T_locked_IBI_tertiles{1, 2}(:, t, ibi, :, 1:300), 1));
    older_ibi_mean(ibi, :, :) = squeeze(mean(HEP_T_locked_IBI_tertiles{2, 2}(:, t, ibi, :, 1:300), 1));

    young_ibi_se(ibi, :, :) = squeeze(std(HEP_T_locked_IBI_tertiles{1, 2}(:, t, ibi, :, 1:300), [], 1))/sqrt(size(HEP_T_locked_IBI_tertiles{1, 2}, 1));
    older_ibi_se(ibi, :, :)= squeeze(std(HEP_T_locked_IBI_tertiles{2, 2}(:, t, ibi, :, 1:300), [], 1))/sqrt(size(HEP_T_locked_IBI_tertiles{1, 2}, 1));
end

    % channels_num = [7, 16, 25, 34];
    channels_num = 51%34;%[7, 25, 43, 57, 60];
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
            plot(x_axis(1:124), squeeze(young_ibi_mean(ibi, p, 1:124)), 'color', [.75 .75 .75],  'LineWidth',2.5)
            hold on
        end
        plot(x_axis, zeros(length(x_axis), 1), '--k', 'LineWidth', 1.5)
%         if strcmp(name_channels, 'EKG')
%             mask = [];
%             if t == 2
%                 mask = find(mask_ecg_simpleRT.mask == 1);
%             elseif t == 3
%                 mask = find(mask_ecg_simpleRT.mask == 1);
%             end
%                 plot(x_axis(mask+125), ones(length(mask), 1)  * -50, 's', 'MarkerSize', 4.5, 'color', 'k', ...
%             'MarkerFaceColor', 'k')
%         elseif t == 2 % simple RT
%             mask = find(mask_simpleRT.mask(p, :) == 1);
%              plot(x_axis(mask+125), ones(length(mask), 1) * -1, 's', 'MarkerSize', 4.5, 'color', 'k', ...
%             'MarkerFaceColor', 'k')
%         elseif t == 3 % gng
%             mask = find(mask_gng.mask(p, :) == 1);
%              plot(x_axis(mask+125), ones(length(mask), 1) * -1, 's', 'MarkerSize', 4.5, 'color', 'k', ...
%             'MarkerFaceColor', 'k')
%         end
        hold off
        box off
        ax = gca;
        if strcmp(name_channels, 'EKG')
            axis([-inf inf -100 250])%axis([-inf .7 -200 200])
        else
            axis([-inf inf -1.50 1.50]);%axis([-inf .7 -1.50 1.50])
        end
%         if t == 1
%             axis([-inf .7 -1.5 1])
%         else
%             axis([-inf .7 -1.5 1])
%         end
        ax.LineWidth = 2.5;
        ax.FontSize = 28;
        ax.FontName = 'Arial';
        xlabel('Time (ms)', 'FontSize', 32, 'FontWeight','normal')
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
            plot(x_axis(1:124), squeeze(older_ibi_mean(ibi, p, 1:124)), 'color', [.75 .75 .75],  'LineWidth',2.5)
            hold on
        end
        plot(x_axis, zeros(length(x_axis), 1), '--k', 'LineWidth', 1.5)
%         if strcmp(name_channels, 'EKG')
%             mask = [];
%             if t == 2
%                 mask = find(mask_ecg_simpleRT.mask == 1);
%             elseif t == 3
%                 mask = find(mask_ecg_simpleRT.mask == 1);
%             end
%                 plot(x_axis(mask+125), ones(length(mask), 1)  * -50, 's', 'MarkerSize', 4.5, 'color', 'k', ...
%             'MarkerFaceColor', 'k')
%         elseif t == 2 % simple RT
%             mask = find(mask_simpleRT.mask(p, :) == 1);
%             plot(x_axis(mask+125), ones(length(mask), 1) * -1, 's', 'MarkerSize', 4.5, 'color', 'k', ...
%             'MarkerFaceColor', [.5 .5 .5])
%         elseif t == 3 % gng
%             mask = find(mask_gng.mask(p, :) == 1);
%              plot(x_axis(mask+125), ones(length(mask), 1) * -1, 's', 'MarkerSize', 4.5, 'color', 'k', ...
%             'MarkerFaceColor', 'k')
%         end
        hold off
        box off
        ax = gca;
        if p == 60 
            axis([-inf  inf -100 250])%axis([-inf .7 -200 200])
        else
            axis([-inf inf  -1.50 1.50])
        end
%         if t == 1
%             axis([-inf inf -1.5 1.5])
%         else
%             axis([-inf inf -2.5 1.5])
%         end
        ax.LineWidth = 2.5;
        ax.FontSize = 28;
        ax.FontName = 'Arial';
        xlabel('Time (ms)', 'FontSize', 32, 'FontWeight','normal')
        ylabel('Amplitude (\muV)', 'FontSize', 32, 'FontWeight','normal')
        name_channels =  channels.chanlocs_with_ekg(p).labels;
        title(['Older - ', tasks{t}, ' - ', name_channels] , 'FontSize', 28, 'FontWeight','normal')

    end

end

%% both groups together
channels = load('chanlocs_with_ekg.mat');
tasks = {'passive', 'simple RT', 'go/no-go'};
load HEP_T_locked_IBI_tertiles;

x_axis=-199:2:400; % in miliseconds

% load time windows where HEP correlates with IBI
mask_dir = 'E:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\LIMO_stats\HEP_TLockedvsIBI\level2_onesample_IBIvsHEP';
mask_simpleRT = load([mask_dir, '\simpleRT\mask']);
mask_gng = load([mask_dir, '\gng\mask']);

for t = 1:length(tasks)

for ibi = 1:3
    ibi_mean(ibi, :, :) = squeeze(mean(cat(1, HEP_T_locked_IBI_tertiles{1, 2}(:, t, ibi, :, 1:300),...
        HEP_T_locked_IBI_tertiles{2, 2}(:, t, ibi, :, 1:300)), 1));

    ibi_se(ibi, :, :) = squeeze(std(cat(1, HEP_T_locked_IBI_tertiles{1, 2}(:, t, ibi, :, 1:300),...
        HEP_T_locked_IBI_tertiles{2, 2}(:, t, ibi, :, 1:300)), [], 1))/sqrt(size(HEP_T_locked_IBI_tertiles{1, 2}, 1));
end

    % channels_num = [7, 16, 25, 34];
    channels_num = 34;%[7, 25, 43, 57, 60];
    colormap cool;
    cmap = colormap;
    colors_ibi = [cmap(1,:); cmap(size(cmap, 1)/2,:); cmap(end,:)];
    for p = channels_num
        name_channels = channels.chanlocs_with_ekg(p).labels;

        figure;
        for ibi = 1:3
            plot(x_axis, squeeze(ibi_mean(ibi, p, :)), 'color', colors_ibi(ibi, :),  'LineWidth',2.5)
            hold on
            jbfill(x_axis, squeeze([ibi_mean(ibi, p, :) + ibi_se(ibi, p, :)])', squeeze([ibi_mean(ibi, p, :) - ibi_se(ibi, p, :)])', colors_ibi(ibi, :), colors_ibi(ibi, :), 1, .2)
            hold on
%             plot(x_axis(1:124), squeeze(ibi_mean(ibi, p, 1:124)), 'color', [.75 .75 .75],  'LineWidth',2.5)
%             hold on
        end
        plot(x_axis, zeros(length(x_axis), 1), '--k', 'LineWidth', 1.5)
        hold off
        box off
        ax = gca;
        if strcmp(name_channels, 'EKG')
            axis([-inf inf -100 250])%axis([-inf .7 -200 200])
        else
            axis([-inf inf -1.50 1.50]);%axis([-inf .7 -1.50 1.50])
        end
        ax.LineWidth = 2.5;
        ax.FontSize = 28;
        ax.FontName = 'Open Sans';
        xlabel('Time (ms)', 'FontSize', 32, 'FontWeight','normal')
        ylabel('Amplitude (\muV)', 'FontSize', 32, 'FontWeight','normal')
        name_channels =  channels.chanlocs_with_ekg(p).labels;
        title([tasks{t}, ' - ', name_channels] , 'FontSize', 28, 'FontWeight','normal')
%         basefilename2 = sprintf('Tepochplot_gng-sRT_young_cut_%d.fig', p);
%         savefig(basefilename2);

    end

end



%% plot effect of IBI on ECG - average of both groups, average of both tasks

load HEP_T_locked_IBI_tertiles;



% load time windows where HEP correlates with IBI
% masks were calculated from 100-700 ms after R peak
mask_ecg_dir = 'E:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\LIMO_stats\ECGvsIBI\level2_onesample_IBIvsECG';
mask_ecg_simpleRT = load([mask_ecg_dir, '\simpleRT\mask']);
mask_ecg_gng = load([mask_ecg_dir, '\gng\mask']);
%     HEP_T_locked_IBI_tertiles = participants x tasks x ibi x channels x time
%avg data across tasks
young_data = squeeze(mean(HEP_T_locked_IBI_tertiles{1, 2}(:, 2:3, :, 60, :), 2));
older_data = squeeze(mean(HEP_T_locked_IBI_tertiles{2, 2}(:, 2:3, :, 60, :), 2));

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
            'MarkerEdgeColor','k','MarkerSize',8, 'LineWidth', 1);
        hold on;
    end

    % axes('XColor','none');
    hold off;
    axis([0 6 -inf inf]);
    xticks([1.5 4.5])
    
    ax = gca;
    ax.XTickLabel= x_label;
    ax.LineWidth = 2.5; 
    ax.FontSize = 24;
    ax.FontName = 'Arial';
    ax.Color = 'none';
%     ax.XTickLabel= [];
    ax.XAxis.FontSize = 24;
    ylabel(y_label_text, 'FontSize', 26, 'FontWeight','normal')
    title(title_text, 'FontSize', 28, 'FontWeight','normal')

%     x0=10;
%     y0=10;
%     width=400;
%     height=400;
%     set(gcf,'position',[x0,y0,width,height])
end