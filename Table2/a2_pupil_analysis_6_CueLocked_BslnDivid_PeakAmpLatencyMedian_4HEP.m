% creating target aligned epochs
% Maria Ribeiro - October 2016 - updated March 2017
clear; close all
younger=[4	9	10	13	15	16	25	26	28	31	33	34	36	42	44	45	46 50 51 53 54 56	59	62 66	68	72	74	76  78  80 81   82  84  85];
older=[7  8	11	12	14	17	19	20	21	22	23	32	35	37	41	43	47 48 49 52 55	57	58	60	61	63	64  65	67	69	70	71	73	75  77 79   83  86];

participants=[younger older];
x=0; y=0;
latency_gng_young=[]; latency_gng_older=[]; latency_simpleRT_young=[]; latency_simpleRT_older=[];
for p=participants

    clearvars -except p x y w z younger older participants...
        BlnDivided_Part_PeakAmpLat_Median_D_CueLocked_YoungGrp BlnDivided_Part_PeakAmpLat_Median_D_CueLocked_OlderGrp...
        BlnDivided_Part_PeakAmpLat_Median_G_CueLocked_YoungGrp BlnDivided_Part_PeakAmpLat_Median_G_CueLocked_OlderGrp...
        latency_gng_young latency_gng_older latency_simpleRT_young latency_simpleRT_older
    
    participant_dir=strcat('M:\ProjectAgingNeuromodulation\AuditoryResearch\AuditoryTask_EyeTracking_EEG_lab94\AB', num2str(p), '\ET');
%     load_directory=strcat(participant_dir, '\ResponseLockedAnalysis');


load([participant_dir filesep 'BaselineDivided_Detection_CueLocked_PD_Correct']);

% load RT and cue-target intervals
load([participant_dir filesep 'D_merged_behavior' filesep 'detection_RT.mat']);
load([participant_dir filesep 'D_merged_behavior' filesep 'warning2imperativetime_detection.mat']);
 
    % find peak amplitude and latency of single trial pupil dilation responses
% for each subject
clear border_trials1 border_trials2 border_trials
[amplitude, latency]=max(BaselineDivided_Detection_CueLocked_PD_Correct(:, 241:1680)');
%delete trials where max was found at border of search window
border_trials1=find(latency==1);
border_trials2=find(latency==1440);%latency set at last point
border_trials=[border_trials1'; border_trials2'];
% tansform latency points in ms
latency_ms=(latency)*1000/240; % zero at target
latency_ms(border_trials)=[];
amplitude(border_trials)=[];

latency_ms=latency_ms';
amplitude=amplitude';

    if isempty(setdiff(p, younger))
        x=x+1;
        BlnDivided_Part_PeakAmpLat_Median_D_CueLocked_YoungGrp(x, :)=[p median(amplitude) median(latency_ms)];
        latency_simpleRT_young=[latency_simpleRT_young; latency_ms];
    elseif isempty(setdiff(p, older))
        y=y+1;
        BlnDivided_Part_PeakAmpLat_Median_D_CueLocked_OlderGrp(y, :)=[p median(amplitude) median(latency_ms)];
        latency_simpleRT_older=[latency_simpleRT_older; latency_ms];
    end


    
    % gng task

load([participant_dir, '\BaselineDivided_GNG_CueLocked_PD_Correct.mat']);

% load RT and cue-target intervals
load([participant_dir filesep 'G_merged_behavior' filesep 'go_RT.mat']);
load([participant_dir filesep 'G_merged_behavior' filesep 'warning2imperativetime_gng.mat']);

        % find peak amplitude and latency of single trial pupil dilation responses
% for each subject
clear border_trials1 border_trials2 border_trials amplitude latency latency_ms
[amplitude, latency]=max(BaselineDivided_GNG_CueLocked_PD_Correct(:, 241:1680)');
%delete trials where max was found at border of search window
border_trials1=find(latency==1);
border_trials2=find(latency==1440);%latency set at last point
border_trials=[border_trials1'; border_trials2'];
% tansform latency points in ms
latency_ms=(latency)*1000/240; % zero at button press
latency_ms(border_trials)=[];
amplitude(border_trials)=[];

latency_ms=latency_ms';
amplitude=amplitude';

    if isempty(setdiff(p, younger))
        BlnDivided_Part_PeakAmpLat_Median_G_CueLocked_YoungGrp(x, :)=[p median(amplitude) median(latency_ms)];
        latency_gng_young=[latency_gng_young; latency_ms];
    elseif isempty(setdiff(p, older))
        BlnDivided_Part_PeakAmpLat_Median_G_CueLocked_OlderGrp(y, :)=[p median(amplitude) median(latency_ms)];
        latency_gng_older=[latency_gng_older; latency_ms];
    end
end


cd('M:\ProjectAgingNeuromodulation\AuditoryResearch\PupilDilation_analysis\BaselineDivided\CueLockedData');
save BlnDivided_Part_PeakAmpLat_Median_G_CueLocked_YoungGrp BlnDivided_Part_PeakAmpLat_Median_G_CueLocked_YoungGrp
save BlnDivided_Part_PeakAmpLat_Median_G_CueLocked_OlderGrp BlnDivided_Part_PeakAmpLat_Median_G_CueLocked_OlderGrp
save BlnDivided_Part_PeakAmpLat_Median_D_CueLocked_YoungGrp BlnDivided_Part_PeakAmpLat_Median_D_CueLocked_YoungGrp
save BlnDivided_Part_PeakAmpLat_Median_D_CueLocked_OlderGrp BlnDivided_Part_PeakAmpLat_Median_D_CueLocked_OlderGrp

figure; 
subplot(2, 2, 1);
hist(latency_simpleRT_young, 20)
ax = gca;
c = ax.Color;
% legend('Detection', 'GNG')
ax.FontSize = 12;
ax.FontName = 'Arial';
ax.Color = 'none';
xlabel('Latency (ms)', 'FontSize', 14, 'FontWeight','bold')
title('simple RT - young', 'FontSize', 14, 'FontWeight','bold')

subplot(2, 2, 2); hist(latency_simpleRT_older, 20)
ax = gca;
c = ax.Color;
% legend('Detection', 'GNG')
ax.FontSize = 12;
ax.FontName = 'Arial';
ax.Color = 'none';
xlabel('Latency (ms)', 'FontSize', 14, 'FontWeight','bold')
title('simple RT - older', 'FontSize', 14, 'FontWeight','bold')

subplot(2, 2, 3); hist(latency_gng_young, 20)
ax = gca;
c = ax.Color;
% legend('Detection', 'GNG')
ax.FontSize = 12;
ax.FontName = 'Arial';
ax.Color = 'none';
xlabel('Latency (ms)', 'FontSize', 14, 'FontWeight','bold')
title('gng - young', 'FontSize', 14, 'FontWeight','bold')

subplot(2, 2, 4); hist(latency_gng_older, 20)
ax = gca;
c = ax.Color;
% legend('Detection', 'GNG')
ax.FontSize = 12;
ax.FontName = 'Arial';
ax.Color = 'none';
xlabel('Latency (ms)', 'FontSize', 14, 'FontWeight','bold')
title('gng - older', 'FontSize', 14, 'FontWeight','bold')

%% plots

figure; 

histogram(latency_simpleRT_young, 40)
hold on
histogram(latency_simpleRT_older, 40)
ax = gca;
c = ax.Color;
% legend('Detection', 'GNG')
ax.FontSize = 12;
ax.FontName = 'Arial';
ax.Color = 'none';
xlabel('Latency (ms)', 'FontSize', 14, 'FontWeight','bold')
title('simple RT', 'FontSize', 14, 'FontWeight','bold')

%% plots

figure; 

histogram(latency_gng_young, 40)
hold on
histogram(latency_gng_older, 40)
ax = gca;
c = ax.Color;
% legend('Detection', 'GNG')
ax.FontSize = 12;
ax.FontName = 'Arial';
ax.Color = 'none';
xlabel('Latency (ms)', 'FontSize', 14, 'FontWeight','bold')
title('gng', 'FontSize', 14, 'FontWeight','bold')