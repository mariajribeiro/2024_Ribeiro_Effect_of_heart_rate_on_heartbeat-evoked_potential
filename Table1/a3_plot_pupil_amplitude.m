% plot amplitude of pupil phasic response
clear;
% eeglab; limo_eeg;

% addpath('C:\eeglab2022.0\plugins\limo_tools-master_Nov2022') 
% check if betas from regression of IBI vs HEP are correlated with
% amplitude of phasic pupil responses
older = [11 12 14 20 21 22 23 32 37 38 41 43 47 48 49 52 55 57 58 63 64 65 67 69 7 70 71 75 8 83 86];
young = [13 15 16 25 26 28 31 33 34 36 4 42 44 45 46 50 51 53 54 56 59 6 62 66 68 72 74 76 78 80 82 84 85 9];
tasks = {'simpleRT', 'gng'};
folder_list = dir(pwd);
load('E:\ProjectAgeingAuditoryTask\heartbeat_evoked_potentials_study\LIMO_stats\expected_chanlocs_both.mat');

% load pupil data
load([pwd, '\pupil_dilation_amplitude\BlnDivided_Part_PeakAmpLat_Median_D_CueLocked_OlderGrp.mat']);
load([pwd, '\pupil_dilation_amplitude\BlnDivided_Part_PeakAmpLat_Median_G_CueLocked_OlderGrp']);
load([pwd, '\pupil_dilation_amplitude\BlnDivided_Part_PeakAmpLat_Median_D_CueLocked_YoungGrp']);
load([pwd, '\pupil_dilation_amplitude\BlnDivided_Part_PeakAmpLat_Median_G_CueLocked_YoungGrp']);

pupil_amplitude{1} = [BlnDivided_Part_PeakAmpLat_Median_D_CueLocked_YoungGrp(:, 1:2); ...
    BlnDivided_Part_PeakAmpLat_Median_D_CueLocked_OlderGrp(:, 1:2)];

pupil_amplitude{2} = [BlnDivided_Part_PeakAmpLat_Median_G_CueLocked_YoungGrp(:, 1:2); ...
    BlnDivided_Part_PeakAmpLat_Median_G_CueLocked_OlderGrp(:, 1:2)];


% limo_random_robust(4,y,X,parameter number,LIMO)
%                    4 = regression analysis
%                    y = data (dim channels, time or freq, subjects)
%                      = data (dim channels, freq, time, subjects)
%                    X = continuous regressor(s)
%                    parameter number = describe which parameter is currently analysed (e.g. 1 - use for maming only)

pupil_young = cell(2, 1); pupil_older = cell(2, 1); 

for t = 1:length(tasks)

    for f = 1:length(folder_list)
        if contains(folder_list(f).name, 'AB')
            % check if participant has pupil data
            subj_idx = find(pupil_amplitude{t}(:, 1) == str2num(folder_list(f).name(3:end)));
            if ~isempty(subj_idx)
                if ismember(str2num(folder_list(f).name(3:end)), young )
                    pupil_young{t} = [pupil_young{t}; str2num(folder_list(f).name(3:end)), pupil_amplitude{t}(subj_idx, 2)];
                else
                    pupil_older{t} = [pupil_older{t}; str2num(folder_list(f).name(3:end)), pupil_amplitude{t}(subj_idx, 2)];
                end
            end
        end
    end
end

%% save data in excel file to perform stats in SPSS
% create table
% participant, group, IBI_avg_passive, IBI_avg_simpleRT, IBI_avg_gng, IBI_std_passive, IBI_std_simpleRT, IBI_std_gng
% IBI_avg_std{group, 2}(old, t, :) = [mean(IBI), std(IBI)];
pupil_table = table([pupil_young{1}(:, 1); pupil_older{1}(:, 1)],...
    [ones(length(pupil_young{1}(:, 1)), 1); ones(length(pupil_older{1}(:, 1)), 1)*2], ...
    [pupil_young{1}(:, 2); pupil_older{1}(:, 2)], ...
    [pupil_young{2}(:, 2); pupil_older{2}(:, 2)], ...
    'VariableNames',{'participant', 'group', 'pupil_amp_simpleRT', 'pupil_amp_gng'});

filename = 'pupil_peak_amplitude.xlsx';
writetable(pupil_table,filename)