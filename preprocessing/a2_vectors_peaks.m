%% R peak

% load the structure WAVES of each eeg run, for example
load ('WAVES_sub-AB76_eeg_sub-AB76_task-simpleRT_run-1_eeg.mat')

% create a variable for the times of R peaks
Rpeak_time=zeros(length(WAVES.QRS),1);

% for subject number 
% 4,6,7,8,9,11,12,13,14,15,16,22,23,26,28,31,32,33,34,36,41,46,47,48,49
% 51,52,53,54,55,57,58,62,63,64,65,66,68,70,72,75,78,80,83,85,86 
% the Rpeak corresponds to WAVES.QRS(:,2)
Rpeak_time=WAVES.QRS(:,2);

% for subject number
% 20,21,25,37,38,42,44,45,50,56,59,67,69,71,74,76,82,84
% the Rpeak corresponds to WAVES.QRS(:,3)
Rpeak_time=WAVES.QRS(:,3);

% compute the times by multiplying the samples per the inverse of the
% sampling rate
for i=1:length(WAVES.QRS)
    Rpeak_time(i,1)=Rpeak_time(i,1)*0.004;
end

% for subjects 4,32,38,49,84 that have ectopic beats, compute a zscore of the 
% RR interval time and remove the R peaks with large zscore
for m=1:length(Rpeak_time)-1
    RR(m)=Rpeak_time(m+1)-Rpeak_time(m);
end
Rscore=zscore(RR);
for n=1:length(Rscore)
    if abs(Rscore(n))>5
        Rpeak_time(n)=0;
        Rpeak_time(n+1)=0;
    end
end
for i=1:length(WAVES.QRS)
    if Rpeak_time(i)~= 0
    Rpeak_time(i,1)=Rpeak_time(i,1)*0.004;
    end
end
Rpeak_time(Rpeak_time==0)=[];

% save the variable
save Rpeak_time_sub-AB76_task-simpleRT_run1 Rpeak_time

% for subject 67 remove wrong R peak
for m=1:length(Rpeak_time)-1
    RR(m)=Rpeak_time(m+1)-Rpeak_time(m);
end
Rscore=zscore(RR);
for n=1:length(Rscore)
    if abs(Rscore(n))>3
        Rpeak_time(n)=0;
        Rpeak_time(n+1)=0;
    end
end
for i=1:length(WAVES.QRS)
    if Rpeak_time(i)~= 0
    Rpeak_time(i,1)=Rpeak_time(i,1)*0.004;
    end
end
Rpeak_time(Rpeak_time==0)=[];


% for subject 76 by looking through eyes the ECG waves, delete elements in the 
% variable Rpeak_time which position is defined by pos
load ('WAVES_sub-AB76_eeg_sub-AB76_task-simpleRT_run-2_eeg.mat')
Rpeak_time=zeros(length(WAVES.QRS),1);
Rpeak_time=WAVES.QRS(:,3);
for i=1:length(WAVES.QRS)
    Rpeak_time(i,1)=Rpeak_time(i,1)*0.004;
end
pos = [1,2,5,8,10,11,14,19,20,23,24,25,27,28,31,35,39,43,44,47,48,50,51,58,72,81,85];
Rpeak_time(pos) = [];
save Rpeak_time_sub-AB76_task-gonogo_run1 Rpeak_time
% same procedure for the other trials and the elements to delete are in the following positions:
pos = [1,24,28,41,49,53,57,61,65,66,68,69,71,74,75,78,79,82,83,91,94,96,100]; % for gng run-2
pos = [12,24,82,89,108,121,128,129,138,153,154,158,159,168,178,186,189,192,195,200,202,203,206,214,215,216,224,230,231,232,236,239,240,245,248,249,251,252,254,256,257,261,263,266,271,275,279,286,287] % for passive run-1
pos = [1,2,3,7,8,9,10,12,16,20,27,34,35,37,38,39,46,49,55,56,59,60,61,62,63,64,75,76,77,78,79,80,83,91,92,94,95,96,99,100,102,103,104] % for simpleRT run-1
pos = [3,6,7,11,12,15,16,20,23,25,26,27,30,32,34,37,41,44,46,50,51,52,57,58,62,64,66,69,71,72,73,77,79,83,85,86,91,94,95,97,98,104,105,108,114,116,118,119,123,126] % for simpleRT run-2



%% T peak

load ('check-WAVES_sub-AB58_eeg_sub-AB58_task-simpleRT_run-1_eeg.mat')
Tpeak_time=zeros(length(WAVES.QRS),1);
Tpeak_time=WAVES.T(:,2);
for i=1:length(WAVES.T)
    Tpeak_time(i,1)=Tpeak_time(i,1)*0.004;
end
Tpeak_time(Tpeak_time==0)=[];
save Tpeak_time_sub-AB58_task-simpleRT_run1 Tpeak_time

% steps for subjects with ectopic beats (same defined for R peaks)
load ('WAVES_sub-AB4_eeg_sub-AB4_task-simpleRT_run-1_eeg.mat')
Tpeak_time=zeros(length(WAVES.T),1);
Tpeak_time=WAVES.T(:,2);
for m=1:length(Tpeak_time)-1
    TT(m)=Tpeak_time(m+1)-Tpeak_time(m);
end
Tscore=zscore(TT);
for n=1:length(Tscore)
    if abs(Tscore(n))>5
        Tpeak_time(n)=0;
        Tpeak_time(n+1)=0;
    end
end
for i=1:length(WAVES.T)
    if Tpeak_time(i)~= 0
    Tpeak_time(i,1)=Tpeak_time(i,1)*0.004;
    end
end
Tpeak_time(Tpeak_time==0)=[];
save Tpeak_time_sub-AB4_task-simpleRT_run1 Tpeak_time

% subjects 22,83 have a negative T peak that is wrong. 
% Change the value of T peak considering the T peak
% equal to T onset detected by the toolbox
load ('WAVES_sub-AB83_eeg_sub-AB83_task-simpleRT_run-2_eeg.mat')
Tpeak_time=zeros(length(WAVES.T),1);
for o=1:length(WAVES.T)
    if WAVES.T(o,2)~=0 && abs(ECG(WAVES.T(o,2)))<ECG(WAVES.T(o,1))
        WAVES.T(o,2)=WAVES.T(o,1);
    end
end
Tpeak_time=WAVES.T(:,2);
for i=1:length(WAVES.T)
    Tpeak_time(i,1)=Tpeak_time(i,1)*0.004;
end
Tpeak_time(Tpeak_time==0)=[];
save Tpeak_time_sub-AB83_task-simpleRT_run2 Tpeak_time

% subjects 53,66 have a negative T peak that is wrong. 
% Change the value of T peak considering the T peak
% equal to T offset detected by the toolbox
load ('WAVES_sub-AB66_eeg_sub-AB66_task-simpleRT_run-2_eeg.mat')
Tpeak_time=zeros(length(WAVES.T),1);
for o=1:length(WAVES.T)
    if ECG(WAVES.T(o,2))<ECG(WAVES.T(o,3))
        WAVES.T(o,2)=WAVES.T(o,3);
    end
end
Tpeak_time=WAVES.T(:,2);
for i=1:length(WAVES.T)
    Tpeak_time(i,1)=Tpeak_time(i,1)*0.004;
end
Tpeak_time(Tpeak_time==0)=[];
save Tpeak_time_sub-AB66_task-simpleRT_run2 Tpeak_time

% remove T peak not properly detected in subjects 34,67,72,76
load ('WAVES_sub-AB43_eeg_sub-AB43_task-gonogo_run-2_eeg.mat')
Tpeak_time=zeros(length(WAVES.T),1);
Tpeak_time=WAVES.T(:,2);
for m=1:length(Tpeak_time)-1
    TT(m)=Tpeak_time(m+1)-Tpeak_time(m);
end
Tscore=zscore(TT);
for n=1:length(Tscore)
    if abs(Tscore(n))>2
        Tpeak_time(n)=0;
        Tpeak_time(n+1)=0;
    end
end
for i=1:length(WAVES.T)
    if Tpeak_time(i)~= 0
    Tpeak_time(i,1)=Tpeak_time(i,1)*0.004;
    end
end
Tpeak_time(Tpeak_time==0)=[];
save Tpeak_time_sub-AB43_task-gonogo_run2 Tpeak_time


