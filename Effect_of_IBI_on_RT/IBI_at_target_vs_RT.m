% calculate HEP 2s interval before cue and 2 s interval after cue
% exclude period of auditory ERP
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

tasks = {'simpleRT', 'gonogo'};

IBI_at_target = cell(2,2);
IBI_at_target_RT_corr = cell(2,2);
IBI_at_target_phase_regression = cell(2,2);
IBI_at_target_phase_interact_regres = cell(2,2);

IBI_at_target_RT_corr_lp = cell(2,2);
IBI_at_target_RT_corr_ep = cell(2,2);

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
    
    IBI_at_target{group, 1} = [IBI_at_target{group, 1}; subj_number];
    IBI_at_target_RT_corr{group, 1} = [IBI_at_target_RT_corr{group, 1}; subj_number];
    IBI_at_target_phase_regression{group, 1} = [IBI_at_target_phase_regression{group, 1}; subj_number];
    IBI_at_target_phase_interact_regres{group, 1} = [IBI_at_target_phase_interact_regres{group, 1}; subj_number];
    IBI_at_target_RT_corr_lp{group, 1} = [IBI_at_target_RT_corr_lp{group, 1}; subj_number];
    IBI_at_target_RT_corr_ep{group, 1} = [IBI_at_target_RT_corr_ep{group, 1}; subj_number];
       
    for t = 1:length(tasks)
        % clear eeglab
        STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];  

        % load data
        filename = {};
        for f = 1:length(folder_content)
            if contains(folder_content(f).name, subjects(s)) && contains(folder_content(f).name, tasks(t))
                filename = folder_content(f).name;
                EEG = pop_loadset('filename', filename, 'filepath', files_folder );
                EEG = pop_select( EEG, 'nochannel',{'EKG'});
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
        
        EEG = pop_epoch( EEG, {  '2'  }, [-2        1], 'newname', 'target_epoch', 'epochinfo', 'yes');
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'gui','off'); 
        EEG = pop_select( EEG, 'time',[ 0 1] );
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 3,'gui','off'); 
        
        % measure RT
        RT = [];
        for ep = 1:length(EEG.epoch)
            responses = find(strcmp(EEG.epoch(ep).eventtype, '5'));
            targets = find(strcmp(EEG.epoch(ep).eventtype, '2'));
            if ~isempty(responses) && ~isempty(targets)
                RT(ep) = EEG.epoch(ep).eventlatency{responses(1)}-EEG.epoch(ep).eventlatency{targets(1)};
            else
                RT(ep) = NaN;
            end
        end
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 4,'retrieve', 4,'study',0); 
        EEG = pop_select( EEG, 'time',[-1 1] );
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 3,'gui','off'); 
            
        % calculate heart rate (or IBI) before target
        IBI_at_target_tmp = []; Rpeaks = [];exclude_epoch = [];
        for ep = 1:length(EEG.epoch)
            Rpeaks = find(strcmp(EEG.epoch(ep).eventtype, 'Rpeak'));
            target = find(strcmp(EEG.epoch(ep).eventtype, '2'));
            if length(Rpeaks) >= 2
                Rb4target = find(cell2mat(EEG.epoch(ep).eventlatency(Rpeaks)) < cell2mat(EEG.epoch(ep).eventlatency(target)));
                Raftertarget = find(cell2mat(EEG.epoch(ep).eventlatency(Rpeaks)) > cell2mat(EEG.epoch(ep).eventlatency(target)));
                if ~isempty(Rb4target) && ~isempty(Raftertarget)
                    IBI_at_target_tmp(ep) = EEG.epoch(ep).eventlatency{Rpeaks(Raftertarget(1))}-EEG.epoch(ep).eventlatency{Rpeaks(Rb4target(end))};
                    cardiac_phase(ep) = (cell2mat(EEG.epoch(ep).eventlatency(target))-EEG.epoch(ep).eventlatency{Rpeaks(Rb4target(end))})/IBI_at_target_tmp(ep)*2*pi;
                else
                    IBI_at_target_tmp(ep) = NaN;
                    cardiac_phase(ep) = NaN;
                end
            else
                IBI_at_target_tmp(ep) = NaN;
                cardiac_phase(ep) = NaN;
%                 exclude_epoch = [exclude_epoch, ep];
            end
        end

        %exclude outlier IBIs
        zscore_IBI = (IBI_at_target_tmp - mean(IBI_at_target_tmp, 'omitnan'))/std(IBI_at_target_tmp, 'omitnan');
        IBI_at_target_tmp(abs(zscore_IBI) > 4) = NaN;
        
         % correlate within participant IBI containing target with RT on a trial-by-trial basis        
        R_IBI_RT = NaN; P_IBI_RT = NaN;
        nan_values = union(find(isnan(IBI_at_target_tmp)), find(isnan(RT)));
        to_include = setdiff(1:length(RT), nan_values);
        [R_IBI_RT,P_IBI_RT] = corrcoef(IBI_at_target_tmp(to_include), RT(to_include));  
        
%         figure; plot(IBI_at_target_tmp(to_include), RT(to_include), 'o')
%         title([num2str(subj_number), ' - ', tasks{t}], 'FontSize', 28, 'FontWeight','normal')
        

%         % regression including as independent variables IBI and cardiac
%         % phase of target onset
%         % [b,bint,r,rint,stats] = regress(y,X)
%         y = zscore(RT(to_include))';
%         x = [ones(length(RT(to_include)), 1), zscore(IBI_at_target_tmp(to_include)'), zscore([sin(cardiac_phase(to_include))]'), zscore([cos(cardiac_phase(to_include))]')];
%         [b,bint,r,rint,stats] = regress(y, x);
        
        
        
%         % regression including as independent variables IBI and cardiac
%         % phase of target onset and interaction terms
%         % [b,bint,r,rint,stats] = regress(y,X)
%         y = zscore(RT(to_include))';
%         x = [ones(length(RT(to_include)), 1), zscore(IBI_at_target_tmp(to_include)'), zscore([sin(cardiac_phase(to_include))]'), zscore([cos(cardiac_phase(to_include))]'), ...
%             zscore(IBI_at_target_tmp(to_include)').*zscore([sin(cardiac_phase(to_include))]'), zscore(IBI_at_target_tmp(to_include)').*zscore([cos(cardiac_phase(to_include))]')];
%         [b_inter,bint,r,rint,stats_inter] = regress(y, x);
        
        
        % correlation between RT and IBI - only trials where target was presented during the
        % second half of cardiac cycle or during the first half of the
        % cardiac cycle
        RT_inc = RT(to_include);
        cardiac_phase_inc = cardiac_phase(to_include);
        IBI_at_target_inc = IBI_at_target_tmp(to_include);
        
        late_phase = find(cardiac_phase_inc > pi);
        [R_IBI_RT_lp,P_IBI_RT_lp] = corrcoef(IBI_at_target_inc(late_phase), RT_inc(late_phase));  
        
        early_phase = find(cardiac_phase_inc < pi);
        [R_IBI_RT_ep,P_IBI_RT_ep] = corrcoef(IBI_at_target_inc(early_phase), RT_inc(early_phase));  
        
%         figure; plot(IBI_at_target_inc(early_phase), RT_inc(early_phase), 'o')
        
        if group == 1
%             IBI_at_target{group, 2}(yng, t) = mean(IBI_at_target_tmp, 'omitnan');
%             IBI_at_target_RT_corr{group, 2}(yng, t, :) = [R_IBI_RT(1, 2),P_IBI_RT(1, 2)];
%             IBI_at_target_phase_regression{group, 2}(yng, t, :) = [b', stats];
%             IBI_at_target_phase_interact_regres{group, 2}(yng, t, :) = [b_inter', stats_inter];
            IBI_at_target_RT_corr_lp{group, 2}(yng, t, :) = [R_IBI_RT_lp(1, 2),P_IBI_RT_lp(1, 2)];
            IBI_at_target_RT_corr_ep{group, 2}(yng, t, :) = [R_IBI_RT_ep(1, 2),P_IBI_RT_ep(1, 2)];
        else
%             IBI_at_target{group, 2}(old, t) = mean(IBI_at_target_tmp, 'omitnan');
%             IBI_at_target_RT_corr{group, 2}(old, t, :) = [R_IBI_RT(1, 2),P_IBI_RT(1, 2)];
%             IBI_at_target_phase_regression{group, 2}(old, t, :) = [b', stats];
%             IBI_at_target_phase_interact_regres{group, 2}(old, t, :) = [b_inter', stats_inter];
            IBI_at_target_RT_corr_lp{group, 2}(old, t, :) = [R_IBI_RT_lp(1, 2),P_IBI_RT_lp(1, 2)];
            IBI_at_target_RT_corr_ep{group, 2}(old, t, :) = [R_IBI_RT_ep(1, 2),P_IBI_RT_ep(1, 2)];
        end
    end       
end

% save IBI_at_target IBI_at_target
% save IBI_at_target_RT_corr IBI_at_target_RT_corr
% save IBI_at_target_phase_regression IBI_at_target_phase_regression

% save IBI_at_target_phase_interact_regres IBI_at_target_phase_interact_regres

save IBI_at_target_RT_corr_lp IBI_at_target_RT_corr_lp
save IBI_at_target_RT_corr_ep IBI_at_target_RT_corr_ep

% save HEP_aftercue HEP_aftercue
% save HEP_b4cue HEP_b4cue
% save HEP_all HEP_all
% save IBI_aftercue IBI_aftercue
% save IBI_b4cue IBI_b4cue


%% plot figures
% 'IBI before and after the cue'
 load IBI_at_target_RT_corr

% plot_2conditions(grp1_cond1, grp1_cond2, grp2_cond1, grp2_cond2, y_label_text, x_label, title_text)
plot_2conditions(IBI_at_target_RT_corr{1, 2}(:,1, 1), IBI_at_target_RT_corr{1, 2}(:,2, 1),...
    IBI_at_target_RT_corr{2, 2}(:,1, 1), IBI_at_target_RT_corr{2, 2}(:,2, 1), 'Correlation coef (r)', {'Simple RT' 'Go/no-go' 'Simple RT' 'Go/no-go'}, '')
hold on
plot([0, 5],zeros(2, 1), '--k')


[h,p,ci,stats] = ttest(IBI_at_target_RT_corr{1, 2}(:,1, 1))
[h,p,ci,stats] = ttest(IBI_at_target_RT_corr{1, 2}(:,2, 1))
[h,p,ci,stats] = ttest(IBI_at_target_RT_corr{2, 2}(:,1, 1))
[h,p,ci,stats] = ttest(IBI_at_target_RT_corr{2, 2}(:,2, 1))

mean(IBI_at_target_RT_corr{1, 2}(:,1, 1))
std(IBI_at_target_RT_corr{1, 2}(:,1, 1))

mean(IBI_at_target_RT_corr{1, 2}(:,2, 1))
std(IBI_at_target_RT_corr{1, 2}(:,2, 1))

mean(IBI_at_target_RT_corr{2, 2}(:,1, 1))
std(IBI_at_target_RT_corr{2, 2}(:,1, 1))

mean(IBI_at_target_RT_corr{2, 2}(:,2, 1))
std(IBI_at_target_RT_corr{2, 2}(:,2, 1))


%% plot correlation coefficients between RT and IBI using only trials where target was presented in the second half of cardiac cycle
 load IBI_at_target_RT_corr_lp

% plot_2conditions(grp1_cond1, grp1_cond2, grp2_cond1, grp2_cond2, y_label_text, x_label, title_text)
plot_2conditions(IBI_at_target_RT_corr_lp{1, 2}(:,1, 1), IBI_at_target_RT_corr_lp{1, 2}(:,2, 1),...
    IBI_at_target_RT_corr_lp{2, 2}(:,1, 1), IBI_at_target_RT_corr_lp{2, 2}(:,2, 1), 'Correlation coef (r)', {'Simple RT' 'Go/no-go' 'Simple RT' 'Go/no-go'}, 'late cardiac phase')
hold on
plot([0, 5],zeros(2, 1), '--k')


[h,p] = ttest(IBI_at_target_RT_corr_lp{1, 2}(:,1, 1))
[h,p] = ttest(IBI_at_target_RT_corr_lp{1, 2}(:,2, 1))
[h,p] = ttest(IBI_at_target_RT_corr_lp{2, 2}(:,1, 1))
[h,p] = ttest(IBI_at_target_RT_corr_lp{2, 2}(:,2, 1))

% both groups together
[h,p] = ttest([IBI_at_target_RT_corr_lp{1, 2}(:,1, 1); IBI_at_target_RT_corr_lp{2, 2}(:,1, 1)])
[h,p] = ttest([IBI_at_target_RT_corr_lp{1, 2}(:,2, 1); IBI_at_target_RT_corr_lp{2, 2}(:,2, 1)])

%% plot correlation coefficients between RT and IBI using only trials where target was presented in the first half of cardiac cycle
 load IBI_at_target_RT_corr_ep

% plot_2conditions(grp1_cond1, grp1_cond2, grp2_cond1, grp2_cond2, y_label_text, x_label, title_text)
plot_2conditions(IBI_at_target_RT_corr_ep{1, 2}(:,1, 1), IBI_at_target_RT_corr_ep{1, 2}(:,2, 1),...
    IBI_at_target_RT_corr_ep{2, 2}(:,1, 1), IBI_at_target_RT_corr_ep{2, 2}(:,2, 1), 'Correlation coef (r)', {'Simple RT' 'Go/no-go' 'Simple RT' 'Go/no-go'}, 'early cardiac phase')
hold on
plot([0, 5],zeros(2, 1), '--k')


[h,p] = ttest(IBI_at_target_RT_corr_ep{1, 2}(:,1, 1))
[h,p] = ttest(IBI_at_target_RT_corr_ep{1, 2}(:,2, 1))
[h,p] = ttest(IBI_at_target_RT_corr_ep{2, 2}(:,1, 1))
[h,p] = ttest(IBI_at_target_RT_corr_ep{2, 2}(:,2, 1))
% both groups together
[h,p] = ttest([IBI_at_target_RT_corr_ep{1, 2}(:,1, 1); IBI_at_target_RT_corr_ep{2, 2}(:,1, 1)])
[h,p] = ttest([IBI_at_target_RT_corr_ep{1, 2}(:,2, 1); IBI_at_target_RT_corr_ep{2, 2}(:,2, 1)])

%%
% 'IBI before and after the cue' - regression coefficient correcting for
% cardiac phase with interaction terms
load IBI_at_target_phase_regression


title_txt = {'constant' 'IBI' 'sin(cardiac phase)' 'cos(cardiac phase)'};

for prd = 1:4
    % plot_2conditions(grp1_cond1, grp1_cond2, grp2_cond1, grp2_cond2, y_label_text, x_label, title_text)
    plot_2conditions(IBI_at_target_phase_regression{1, 2}(:,1, prd), IBI_at_target_phase_regression{1, 2}(:,2, prd),...
        IBI_at_target_phase_regression{2, 2}(:,1, prd), IBI_at_target_phase_regression{2, 2}(:,2, prd), 'Regression coefficient', {'Simple RT' 'Go/no-go' 'Simple RT' 'Go/no-go'}, title_txt{prd})
    hold on
    plot([0, 5],zeros(2, 1), '--k')
end

[h,p] = ttest(IBI_at_target_phase_regression{1, 2}(:,1, 2))
[h,p] = ttest(IBI_at_target_phase_regression{1, 2}(:,2, 2))
[h,p] = ttest(IBI_at_target_phase_regression{2, 2}(:,1, 2))
[h,p] = ttest(IBI_at_target_phase_regression{2, 2}(:,2, 2))


[h,p] = ttest(IBI_at_target_phase_regression{1, 2}(:,1, 3))
[h,p] = ttest(IBI_at_target_phase_regression{1, 2}(:,2, 3))
[h,p] = ttest(IBI_at_target_phase_regression{2, 2}(:,1, 3))
[h,p] = ttest(IBI_at_target_phase_regression{2, 2}(:,2, 3))

[h,p] = ttest(IBI_at_target_phase_regression{1, 2}(:,1, 4))
[h,p] = ttest(IBI_at_target_phase_regression{1, 2}(:,2, 4))
[h,p] = ttest(IBI_at_target_phase_regression{2, 2}(:,1,4))
[h,p] = ttest(IBI_at_target_phase_regression{2, 2}(:,2, 4))




%%
% 'IBI before and after the cue' - regression coefficient correcting for
% cardiac phase with interaction terms
load IBI_at_target_phase_interact_regres


title_txt = {'constant' 'IBI' 'sin(cardiac phase)' 'cos(cardiac phase)' 'interc sin' 'interc cos'};

for prd = 1:6
    % plot_2conditions(grp1_cond1, grp1_cond2, grp2_cond1, grp2_cond2, y_label_text, x_label, title_text)
    plot_2conditions(IBI_at_target_phase_interact_regres{1, 2}(:,1, prd), IBI_at_target_phase_interact_regres{1, 2}(:,2, prd),...
        IBI_at_target_phase_interact_regres{2, 2}(:,1, prd), IBI_at_target_phase_interact_regres{2, 2}(:,2, prd), 'Regression coefficient', {'Simple RT' 'Go/no-go' 'Simple RT' 'Go/no-go'}, title_txt{prd})
    hold on
    plot([0, 5],zeros(2, 1), '--k')
end
%%
[h,p] = ttest(IBI_at_target_phase_interact_regres{1, 2}(:,1, 2))
[h,p] = ttest(IBI_at_target_phase_interact_regres{1, 2}(:,2, 2))
[h,p] = ttest(IBI_at_target_phase_interact_regres{2, 2}(:,1, 2))
[h,p] = ttest(IBI_at_target_phase_interact_regres{2, 2}(:,2, 2))


[h,p] = ttest(IBI_at_target_phase_interact_regres{1, 2}(:,1, 3))
[h,p] = ttest(IBI_at_target_phase_interact_regres{1, 2}(:,2, 3))
[h,p] = ttest(IBI_at_target_phase_interact_regres{2, 2}(:,1, 3))
[h,p] = ttest(IBI_at_target_phase_interact_regres{2, 2}(:,2, 3))

[h,p] = ttest(IBI_at_target_phase_interact_regres{1, 2}(:,1, 4))
[h,p] = ttest(IBI_at_target_phase_interact_regres{1, 2}(:,2, 4))
[h,p] = ttest(IBI_at_target_phase_interact_regres{2, 2}(:,1,4))
[h,p] = ttest(IBI_at_target_phase_interact_regres{2, 2}(:,2, 4))

[h,p] = ttest(IBI_at_target_phase_interact_regres{1, 2}(:,1, 5))
[h,p] = ttest(IBI_at_target_phase_interact_regres{1, 2}(:,2, 5))
[h,p] = ttest(IBI_at_target_phase_interact_regres{2, 2}(:,1, 5))
[h,p] = ttest(IBI_at_target_phase_interact_regres{2, 2}(:,2, 5))

[h,p] = ttest([IBI_at_target_phase_interact_regres{1, 2}(:,2, 5); IBI_at_target_phase_interact_regres{2, 2}(:,2, 5)])


[h,p] = ttest(IBI_at_target_phase_interact_regres{2, 2}(:,2, 5))

%% p value of full regression model
% plot_2conditions(grp1_cond1, grp1_cond2, grp2_cond1, grp2_cond2, y_label_text, x_label, title_text)
plot_2conditions(IBI_at_target_phase_interact_regres{1, 2}(:,1,end-1), IBI_at_target_phase_interact_regres{1, 2}(:,2, end-1),...
    IBI_at_target_phase_interact_regres{2, 2}(:,1, end-1), IBI_at_target_phase_interact_regres{2, 2}(:,2, end-1), 'P value', {'Simple RT' 'Go/no-go' 'Simple RT' 'Go/no-go'}, '')
hold on
plot([0, 5],zeros(2, 1), '--k')
%%
plot_2conditions(IBI_at_target_phase_interact_regres{1, 2}(:,1,7), IBI_at_target_phase_interact_regres{1, 2}(:,2, 7),...
    IBI_at_target_phase_interact_regres{2, 2}(:,1, 7), IBI_at_target_phase_interact_regres{2, 2}(:,2, 7), 'R-square', {'Simple RT' 'Go/no-go' 'Simple RT' 'Go/no-go'}, '')
hold on
plot([0, 5],zeros(2, 1), '--k')



%%
plot_2conditions(IBI_at_target_phase_regression{1, 2}(:,1, 3), IBI_at_target_phase_regression{1, 2}(:,2, 3),...
    IBI_at_target_phase_regression{2, 2}(:,1, 3), IBI_at_target_phase_regression{2, 2}(:,2, 3), 'Coefficients', {'Simple RT' 'Go/no-go' 'Simple RT' 'Go/no-go'}, 'sin')
hold on
plot([0, 5],zeros(2, 1), '--k')

[h,p] = ttest(IBI_at_target_phase_regression{1, 2}(:,1, 4))
[h,p] = ttest(IBI_at_target_phase_regression{1, 2}(:,2, 4))
[h,p] = ttest(IBI_at_target_phase_regression{2, 2}(:,1, 4))
[h,p] = ttest(IBI_at_target_phase_regression{2, 2}(:,2, 4))

%%
plot_2conditions(IBI_at_target_phase_regression{1, 2}(:,1, 4), IBI_at_target_phase_regression{1, 2}(:,2, 4),...
    IBI_at_target_phase_regression{2, 2}(:,1, 4), IBI_at_target_phase_regression{2, 2}(:,2, 4), 'Coefficients', {'Simple RT' 'Go/no-go' 'Simple RT' 'Go/no-go'}, 'cos')
hold on
plot([0, 5],zeros(2, 1), '--k')

%% plot_2conditions(grp1_cond1, grp1_cond2, grp2_cond1, grp2_cond2, y_label_text, x_label, title_text)
plot_2conditions(IBI_at_target_phase_regression{1, 2}(:,1, 7), IBI_at_target_phase_regression{1, 2}(:,2, 7),...
    IBI_at_target_phase_regression{2, 2}(:,1, 7), IBI_at_target_phase_regression{2, 2}(:,2, 7), 'p value', {'Simple RT' 'Go/no-go' 'Simple RT' 'Go/no-go'}, '')
hold on
plot([0, 5],zeros(2, 1), '--k')

%% functions
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
            'MarkerEdgeColor', 'k','MarkerSize',8, 'LineWidth', 1);
        hold on;
    end

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

function plot_timecourse(data1_mean, data1_se, data2_mean, data2_se, color1, color2, title_txt)
        figure;
        plot(x_axis, data1_mean, 'color', color1,  'LineWidth',1.5) %plot(light blue)
        hold on
        jbfill(x_axis, [data1_mean + data1_se], [data1_mean - data1_se], color1, color1)
        hold on
        plot(x_axis, data2_mean, '--r',  'LineWidth',1.5) 
        hold on
        jbfill(x_axis, [data2_mean + data2_se], [data2_mean - data2_se], color2, color2)
        hold on
        plot(x_axis, zeros(length(x_axis), 1), ':k')
        hold off
        box off
        ax = gca;
        axis([-.150 .650 -1.5 1.5])%axis([-inf inf -1.5 1.5])
        ax.LineWidth = 2.5;
        ax.FontSize = 28;
        ax.FontName = 'Arial';
        xlabel('Time (s)', 'FontSize', 32, 'FontWeight','normal')
        ylabel('Amplitude (\muV)', 'FontSize', 32, 'FontWeight','normal')
        title(title_txt , 'FontSize', 32, 'FontWeight','normal')
end