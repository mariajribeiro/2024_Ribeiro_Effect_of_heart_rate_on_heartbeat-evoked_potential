% Data from directory
dirData = dir('C:\Users\april\Desktop\erasmus - tesi magistrale\ECG_Toolbox_DEI\ecgSimulator3.0\ecgSimulator3.0\DATA\DATA_eeg');      %# Get the data for the current directory
dirIndex = [dirData.isdir];

% List of the subdirectories
subDirs = {dirData(dirIndex).name};  

% Index of subdirectories that are not '.' or '..'
validIndex = ~ismember(subDirs,{'.','..'});  
                                              
% Loop over valid subdirectories                                               
  for iDir = find(validIndex)
    % Subdirectory path  
    nextDir = fullfile('C:\Users\april\Desktop\erasmus - tesi magistrale\ECG_Toolbox_DEI\ecgSimulator3.0\ecgSimulator3.0\DATA\DATA_eeg',subDirs{iDir});    
    nextdirData = dir(nextDir);
    nextdirIndex = [nextdirData.isdir];
       
    % List of the files
    fileList = {nextdirData(~nextdirIndex).name}';
  
    % Loop over eeg runs for each participant
    for k = 1:length(fileList)
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
    EEG = pop_loadset('filename',fileList(k));
    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
    EEG = eeg_checkset( EEG );
    
    clear start_latency end_latency
    
    % remove data before event 9 if exists
    for evn = 1:length(EEG.event)
        if EEG.event(evn).type == 9
           start_latency = EEG.event(evn).latency;
        elseif EEG.event(evn).type == 10
               end_latency = EEG.event(evn).latency;
        end
    end
    A = exist ('start_latency');
    if A == 0;
       start_latency = 0;
    end
 
    %create epoch between event 9 and event 10
    N_samples = end_latency-start_latency;
    end_epoch = N_samples*0.004;
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = eeg_checkset( EEG );
    if start_latency == 0
       EEG = pop_epoch( EEG, [0      'end_epoch'], 'newname', 'data pruned with ICA epochs', 'epochinfo', 'yes');
       [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'gui','off'); 
       EEG = eeg_checkset( EEG );
    else
       EEG = pop_epoch( EEG, {  '9'  }, [0      'end_epoch'], 'newname', 'data pruned with ICA epochs', 'epochinfo', 'yes');
       [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'gui','off'); 
       EEG = eeg_checkset( EEG );
    end
    
    
    % select EKG and resample 
    EEG = pop_select( EEG, 'channel',{'EKG'});
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'gui','off'); 
    EEG = eeg_checkset( EEG );
    EEG = pop_resample( EEG, 250);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'gui','off');
     ecgdata = EEG.data;

     % Path Initialization
     load OPTION
     global OPTION
     rootCV=[char(pwd) '\'];
     ecgPathConfig(rootCV);
     display  =0;
     maximise =OPTION.M;
     default=0;

     % Load and PreProcessing
     DAT.File = double(ecgdata);
     DAT.Samp = 250;
     DAT.Ini = 1;
     DAT.End = length(ecgdata);
     datFile=DAT.File;
     datSamp=DAT.Samp;
     datIni =DAT.Ini;
     datEnd =DAT.End;
     datChannel = 1; 
     
     % Normalization
     norSamp    = 250;          % mSec   :Normalized Samplig rate
     norWin     = 10;           % Sec    :Window
     norMinV    =-1;           % Num    :Minimum value (after normalization)
     norMaxV    = 1;            % Num    :Maximum value (after normalization)
     norType    = 0;            % 0/1    0: [norMinV .. norMaxV]
     %                                   1: mean(X)=0;  std(X)=1

     % Load data
     DAT.File   = datFile;
     DAT.Channel= datChannel;
     DAT.Samp   = datSamp;
     DAT.Ini    = datIni;
     DAT.End    = datEnd;
     DAT.ECG    = zeros(datEnd-datIni+1,1);

     % Load ECG signal  
     ECG = DAT.File;
     tECG = 0:(1/DAT.Samp):length(ECG)/DAT.Samp;  % Time vector (seconds)
     tECG = tECG(1:end-1);

     % Preprocessing
     disp('~~ BEGIN : Preprocessing ')
     ECG= ecgRemoveBaseLine(ECG, norSamp, display);
     ECG= ecgNoiseReduction(ECG, norSamp, 20, display);  % LowPass 20 Hz
     ECG= ecgNormalize( ECG, norMinV, norMaxV, norType);
     disp('~~ END   : Preprocessing ')
    
    if iDir == 33 || iDir == 42 || iDir == 55 || iDir == 56 || iDir == 58    %subject 51,6,74,75,78
     ECG= -ECG;
    end
     
     % Segmentation
%  Detection of QRS peaks, Q wave onset and S wave offset
%  WAVES= ECG_Segment(Signal, Sample_Freq, Window_Size_Sec, Display)
%  Signal       	:ECG signal vector			[vector of double]
%  Sample_Freq      :Ecg signal sampling frequency	in miliseconds [int]
%  Window_Size_Sec  :Length of the reading window in seconds	  [int]
%  Display       	:0/1 display or not the ecg and segmentation	  [int]	
%¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨
%  WAVES. {structure}   onset, peak and offset of the P, QRS and T waves
%       .P         	:Pwave -  (ponset ppeak poffset) 	[matrix of double]
%       .QRS        :QRS complex – (qonset qpeak rpeak speak soffset)
%											[matrix of double]
%       .T          :Twave – (tonset tpeak toffset) 	[matrix of double]
%                 	  NOTE- Absent waves have zero index
%¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨

    disp('~~ BEGIN : ECG Segmentation ')
    WAVES = ECG_Segment(ECG, norSamp, norWin, display);
    disp('~~ END   : ECG Segmentation ')

    % some ecg have some R peaks that coincides with the R onset 
    % subjects 11,12,13,16,23,36,41,46,47,48,52,54,57,58,63,65,68,70,76,78,86
    if iDir == 3 || iDir == 4 || iDir == 5 || iDir == 8 || iDir == 12 || iDir == 20 || iDir == 24 || iDir == 28 || iDir == 29 || iDir == 30 || iDir == 34 || iDir == 36 || iDir == 39 || iDir == 40 || iDir == 44 || iDir == 46 || iDir == 49 || iDir == 52 || iDir == 57 || iDir == 58 || iDir == 65
    for m = 1:length(WAVES.QRS)
        if ECG(WAVES.QRS(m,3))>ECG(WAVES.QRS(m,2))
           WAVES.QRS(m,2) = WAVES.QRS(m,3);
           WAVES.QRS(m,3) = WAVES.QRS(m,4);
           WAVES.QRS(m,4) = WAVES.QRS(m,5);
        end
    end
    end
    % ecg has some R peaks that coincides with the S onset
    % sub AB14
    if iDir == 6
      for m = 1:length(WAVES.QRS)
          if ECG(WAVES.QRS(m,3))>ECG(WAVES.QRS(m,2))
             WAVES.QRS(m,5) = WAVES.QRS(m,4);
             WAVES.QRS(m,4) = WAVES.QRS(m,3);
             WAVES.QRS(m,3) = WAVES.QRS(m,2);
             WAVES.QRS(m,2) = WAVES.QRS(m,1);
          end
      end
    end
    
    % Save WAVES in the same folder of eegdata
    basefilename1 = sprintf('WAVES_%s.mat', fileList{k});
    fullmatfile1 = fullfile(nextDir, basefilename1);
    save(fullmatfile1,'WAVES');
    
    % Plot WAVES
    WAVES=[ WAVES.P WAVES.QRS WAVES.T];
    
    N=length(ECG);
    t=0:(1/DAT.Samp):N/DAT.Samp;
    t=t(1:end-1);

    figure;
    if maximise
       ecgMaximize;
    end

    %subplot(211)
    lineW=1.0;
    subplot(211)
    plot(ECG,'g');
    ecgLabel('Signal','Index', 'Ecg')
    zoom on
    hold on

    % Rpeak
    plot(WAVES(:,6),ECG(WAVES(:,6)),'+r');

    QwaveB=WAVES(:,4);
    TwaveE=WAVES(:,8);
    for i=1:length(QwaveB)
        plot(QwaveB(i):TwaveE(i), ECG(QwaveB(i):TwaveE(i)),'c', 'Linewidth',lineW);
    end
    plot(WAVES(:,4),ECG(WAVES(:,4)),'.k');
    plot(WAVES(:,5),ECG(WAVES(:,5)),'ob');
    plot(WAVES(:,6),ECG(WAVES(:,6)),'+r');
    plot(WAVES(:,7),ECG(WAVES(:,7)),'ob');
    plot(WAVES(:,8),ECG(WAVES(:,8)),'.k');

    % Twave
    TwaveB=WAVES(WAVES(:,9)~=0,9);
    TwaveE=WAVES(WAVES(:,11)~=0,11);
    for i=1:length(TwaveB)
        plot(TwaveB(i):TwaveE(i), ECG(TwaveB(i):TwaveE(i)),'r', 'Linewidth',lineW);
    end
    plot(WAVES(WAVES(:,9)~=0,9),ECG(WAVES(WAVES(:,9)~=0,9)),'.r');
    plot(WAVES(WAVES(:,10)~=0,10),ECG(WAVES(WAVES(:,10)~=0,10)),'or');
    plot(WAVES(WAVES(:,11)~=0,11),ECG(WAVES(WAVES(:,11)~=0,11)),'.r');

    % Pwave
    PwaveB=WAVES(WAVES(:,1)~=0,1);
    PwaveE=WAVES(WAVES(:,3)~=0,3);
    for i=1:length(PwaveB) 
        plot(PwaveB(i):PwaveE(i), ECG(PwaveB(i):PwaveE(i)),'b', 'Linewidth',lineW);
    end
    plot(WAVES(WAVES(:,1)~=0,1),ECG(WAVES(WAVES(:,1)~=0,1)),'.b');
    plot(WAVES(WAVES(:,2)~=0,2),ECG(WAVES(WAVES(:,2)~=0,2)),'ob');
    plot(WAVES(WAVES(:,3)~=0,3),ECG(WAVES(WAVES(:,3)~=0,3)),'.b');

    %subplot(212)
    lineW=1.5;
    subplot(212)
    hold on;
    plot(ECG,'g')
    ecgLabel('','', 'Ecg')
    zoom on
    
    % Rpeak
    plot(WAVES(:,6),ECG(WAVES(:,6)),'+r');

    QwaveB=WAVES(:,4);
    TwaveE=WAVES(:,8);
    for i=1:length(QwaveB)
        plot(QwaveB(i):TwaveE(i), ECG(QwaveB(i):TwaveE(i)),'c', 'Linewidth',lineW);
    end
    plot(WAVES(:,4),ECG(WAVES(:,4)),'.k');
    plot(WAVES(:,5),ECG(WAVES(:,5)),'ob');
    plot(WAVES(:,6),ECG(WAVES(:,6)),'+r');
    plot(WAVES(:,7),ECG(WAVES(:,7)),'ob');
    plot(WAVES(:,8),ECG(WAVES(:,8)),'.k');
 
    % Twave
    TwaveB=WAVES(WAVES(:,9)~=0,9);
    TwaveE=WAVES(WAVES(:,11)~=0,11);
    for i=1:length(TwaveB)
        plot(TwaveB(i):TwaveE(i), ECG(TwaveB(i):TwaveE(i)),'r', 'Linewidth',lineW);
    end
    plot(WAVES(WAVES(:,9)~=0,9),ECG(WAVES(WAVES(:,9)~=0,9)),'.r');
    plot(WAVES(WAVES(:,10)~=0,10),ECG(WAVES(WAVES(:,10)~=0,10)),'or');
    plot(WAVES(WAVES(:,11)~=0,11),ECG(WAVES(WAVES(:,11)~=0,11)),'.r');

    % Pwave
    PwaveB=WAVES(WAVES(:,1)~=0,1);
    PwaveE=WAVES(WAVES(:,3)~=0,3);
    for i=1:length(PwaveB) 
        plot(PwaveB(i):PwaveE(i), ECG(PwaveB(i):PwaveE(i)),'b', 'Linewidth',lineW);
    end
    plot(WAVES(WAVES(:,1)~=0,1),ECG(WAVES(WAVES(:,1)~=0,1)),'.b');
    plot(WAVES(WAVES(:,2)~=0,2),ECG(WAVES(WAVES(:,2)~=0,2)),'ob');
    plot(WAVES(WAVES(:,3)~=0,3),ECG(WAVES(WAVES(:,3)~=0,3)),'.b');

    % Add slider
    zoom on
    dx=DAT.Samp*5; % 5 seconds
    dx=min(dx,length(ECG));
    a=gca;
    b=gcf;
    set(gcf,'doublebuffer','on');
    set(a,'xlim',[0 dx]);
    set(a,'ylim',[min(ECG) max(ECG) ]);
    pos=get(a,'position');
    Newpos=[pos(1) pos(2)-0.06 pos(3) 0.03];    
    xmax=length(ECG);
    S=['set(gca,''xlim'',get(gcbo,''value'')+[0 ' num2str(dx) '])'];
    hui=uicontrol('style','slider',...
    'units','normalized','position',Newpos,...
    'callback',S,'min',0,'max',xmax-dx);

    % Save plot in the same folder of eegdata
    basefilename2 = sprintf('plot_WAVES_%s.fig', fileList{k});
    fullmatfile2 = fullfile(nextDir, basefilename2);
    savefig(fullmatfile2);
    end
    end
