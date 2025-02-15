function outliers_evnt = find_cardiac_cycle_outliers(EEG)

    % calculate heart rate for each cycle - to exclude cycles with
    % ectopic or wrongly assigned peaks - toolbox mistakes
    IBI_tmp = []; RT_interval= [];
    R_latency = []; T_latency = [];
    for evnt = 1:length(EEG.event)
        if strcmp(EEG.event(evnt).type, 'Rpeak')
            R_latency = [R_latency; EEG.event(evnt).latency*1000/EEG.srate evnt];
        elseif strcmp(EEG.event(evnt).type, 'Tpeak')
            T_latency = [T_latency; EEG.event(evnt).latency*1000/EEG.srate evnt];
        end
    end

    % create matrix with for each line R latency, corresponding event number for R peak, RR interval, 
    % RT interval, and corresponding event number for T peak

    cardiac_measures = [];

    for r = 1:length(R_latency)-1
        cardiac_measures(r, 1) = R_latency(r, 1);
        cardiac_measures(r, 2) = R_latency(r, 2); %event number
        cardiac_measures(r, 3) = R_latency(r+1, 1)-R_latency(r, 1); %IBI
        if ~isempty(T_latency)
            Tpeak = find(T_latency(:, 1) > R_latency(r, 1) & T_latency(:, 1) < R_latency(r+1, 1));
        else
            Tpeak = [];
        end
        if length(Tpeak) == 1
            cardiac_measures(r, 4) = T_latency(Tpeak, 1) - R_latency(r, 1);
            cardiac_measures(r, 5) = T_latency(Tpeak, 2); %event number
        else
            cardiac_measures(r, 4) = NaN;
            cardiac_measures(r, 5) = NaN;
        end

    end


    % find outliers - ectopic beats or cycles where ECG was not
    % properly segmented
    RR_zscore = zscore(cardiac_measures(:, 3));
    RT_zscore = (cardiac_measures(:, 4) - mean(cardiac_measures(:, 4), 'omitnan'))/std(cardiac_measures(:, 4), 'omitnan');

    outliers_RR = find(abs(RR_zscore) > 4);
    outliers_RT = find(abs(RT_zscore) > 4);


    outliers = unique([outliers_RR; outliers_RT; outliers_RR+1]);
    outliers_evnt = [cardiac_measures(outliers, 2); cardiac_measures(outliers, 5)]';

    outliers_evnt(isnan(outliers_evnt)) = [];
    
%     IBI = cardiac_measures(:, 3);
%     IBI(outliers) = [];

%     % delete events associated with outlier cycles
%     EEG = pop_selectevent( EEG, 'omitevent', outliers_evnt ,'deleteevents','on');

end