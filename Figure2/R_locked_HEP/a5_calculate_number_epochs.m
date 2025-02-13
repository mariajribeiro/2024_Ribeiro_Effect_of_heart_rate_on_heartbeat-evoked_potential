% run second level LIMO stats
clear;

older = [11 12 14 20 21 22 23 32 37 38 41 43 47 48 49 52 55 57 58 63 64 65 67 69 7 70 71 75 8 83 86];
young = [13 15 16 25 26 28 31 33 34 36 4 42 44 45 46 50 51 53 54 56 59 6 62 66 68 72 74 76 78 80 82 84 85 9];
% tasks = {'simpleRT', 'gng'}; 
folder_list = dir(pwd);
% 300 ms epochs from 330 up to 630 after R peak - same size epoch as teh
% T-locked analyses

% exclude participants that were not included in T peak analyses: 21 50 55
% - T peak not detected [21, 50, 55]

epoch_RHEP_young = []; epoch_RHEP_older = [];
for f = 1:length(folder_list)
    if contains(folder_list(f).name, 'AB') && ~ismember(str2double(folder_list(f).name(3:end)), [21, 50, 55])
        load([folder_list(f).folder, filesep, folder_list(f).name, filesep, 'categorical_variable.mat'])
        if ismember(str2num(folder_list(f).name(3:end)), young)
            epoch_RHEP_young = [epoch_RHEP_young; length(find(categorical_variable==1)), length(find(categorical_variable==2))];
        elseif ismember(str2num(folder_list(f).name(3:end)), older)
            epoch_RHEP_older = [epoch_RHEP_older; length(find(categorical_variable==1)), length(find(categorical_variable==2))];
        end
    end
end

young_mean = mean(epoch_RHEP_young, 1)
young_std = std(epoch_RHEP_young, [], 1)

older_mean = mean(epoch_RHEP_older, 1)
older_std = std(epoch_RHEP_older, [], 1)