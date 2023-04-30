function lfpShft = ZavSynchLFP(zavp, hd, segms, segmEdge, lfpEx, rCh, rawData, nlxVer)
%lfpShft = ZavSynchLFP(zp, hd, segms, segmEdge, lfpEx, rCh, rawData)
%cut lfp traces phased with respect to stimulus moments
%
%INPUTS
%zp - structure with parameters
%hd - header of data file
%segms - matrix of wanted segments start points, in from:
%        segms(:, 1) - position of synchro-event (ms)
%        segms(:, 2) - number of channel where syncro-event was detected
%        segms(:, 3) - number of sweeps where syncro-event was detected
%        segms(:, 4) - number of stimulus (in matrix zp.realStim) or number of trough in brst matrix
%segmEdge - left and right shifts from synchro-point (ms)
%           full view interval is between segms(sn, 1) + segmEdge(1) and segms(sn, 1) + segmEdge(2)
%lfpEx - external lfp
%rCh - channels to be read
%rawData - boolean variable. If raw data needed rawData = 1(true)
%nlxVer - version ov NLX data (0 - initial, 1 - new, ...)
%
%OUTPUTS
%lfpShft - lfp phased with respect to stimulus moments
%

if rawData %read raw data
    sn = 1;
    pL = segms(sn, 1) + segmEdge(1);%start time (ms from record begin), left point
    pR = segms(sn, 1) + segmEdge(2);%end time (ms from record begin), right point
    
    dpL = 0;%points befor record
    if (pL < 0) %wrong start time
        dpL = round(-pL * 1e3 / hd.ch_si(rCh(1)));%points before record
        pL = 0;%load frome begin (corrected start time, ms)
    end
    recLen = hd.recChDurations(rCh(1)) * 1e-3;%raw signal length (ms)
    dpR = 0;%points after record
    if (pR > recLen) %wrong end time
        dpR = round((pR - recLen) * 1e3 / hd.ch_si(rCh(1)));%points after record
        pR = recLen;%load up to end (corrected stop time, ms)
    end
    
    if ((((pL >= 0) && (pL <= recLen)) || ((pR >= 0) && (pR <= recLen))) && (pL <= pR)) %at least part of data exist
        lfpShft = ZavLoadData(zavp.file, hd, rCh, segms(sn, 3), pL, pR, nlxVer);%load the first sweep from the first channel
        lfpShft = [zeros(dpL, size(lfpShft, 2)); lfpShft; zeros(dpR, size(lfpShft, 2))];%add samples before and after record
    else
        lfpShft = zeros(round((pR - pL) * 1e3 / hd.ch_si(rCh(1))) + dpL + dpR, numel(rCh));%no data from record
    end
    segmLen = size(lfpShft, 1);%length of needed segments (samples)
    lfpShft = cat(3, lfpShft, zeros(segmLen, length(rCh), size(segms, 1) - 1));%lfp phased with respect to stimuli moments
    
    for sn = 2:size(segms, 1) %run over segments
        pL = segms(sn, 1) + segmEdge(1);%start time (ms from record begin)
        pR = segms(sn, 1) + segmEdge(2);%end time (ms from record begin)
        
        dpL = 0;%points befor record
        if (pL < 0) %wrong start time
            dpL = -pL / hd.ch_si(rCh(1));%points before record
            pL = 0;%load frome begin (corrected start time, ms)
        end
        dpR = 0;%points after record
        if (dpR > hd.dataPtsPerChan(rCh(1))) %wrong end time
            dpR = (pR - hd.dataPtsPerChan(rCh(1))) / hd.ch_si(rCh(1));%points after record
            pR = hd.dataPtsPerChan(rCh(1));%load up to end (corrected stop time, ms)
        end
        
        if ((((pL >= 0) && (pL <= recLen)) || ((pR >= 0) && (pR <= recLen))) && (pL <= pR)) %at least part of data exist
            lfpShft(:, :, sn) = [zeros(dpL, size(lfpShft, 2)); ...
                                 ZavLoadData(zavp.file, hd, rCh, segms(sn, 3), pL, pR, nlxVer); ... raw data phased with respect to synchro-event
                                 zeros(dpR, size(lfpShft, 2))];
        else
            lfpShft(:, :, sn) = [zeros(dpL, size(lfpShft, 2)); ...
                                 zeros(round((pR - pL) * 1e3 / hd.ch_si(rCh(1))) + dpL + dpR, numel(rCh)); ... no data from record
                                 zeros(dpR, size(lfpShft, 2))];
        end
    end
else %read resampled data
    segms(:, 1) = zavp.dwnSmplFrq * segms(:, 1) / 1e3;%convert ms to samples
    segmEdge = round(zavp.dwnSmplFrq * segmEdge / 1e3);%convert ms to samples
    segms(:, 1) = round(segms(:, 1)) + 1;%position of synchro-impuls (samples)
    segmLen = diff(segmEdge) + 1;%length of needed segments (samples)
    lfpShft = zeros(segmLen, length(rCh), size(segms, 1));%lfp phased with respect to stimuli moments

    dataLen = size(lfpEx, 1);%length of part of data
    for sn = 1:size(segms, 1) %run over segments
        siPsn = segms(sn, 1);%position of synchro-impuls (samples)
        k = siPsn + segmEdge(1);%negative number equal to number of points to be skipped
        pL = max(1, k);%first point to be read from lfp
        k = ((1 - k) * (k <= 0)) + 1;%first point in lfpShft for to write
        pR = min(siPsn + segmEdge(2), dataLen);%last point to be read from lfp
        t = siPsn + segmEdge(2) - dataLen;%negative number equal to number of points to be skipped
        t = t * (t >= 0);%(size(lfpShft, 1) - t) - last point in lfpShft for write
        lfpShft(k:(end - t), :, sn) = lfpEx(pL:pR, rCh, segms(sn, 3));%resampled data phased with respect to synchro-event
    end
end

