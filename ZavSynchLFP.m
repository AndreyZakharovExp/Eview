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
%           full view interval is between Si + segmEdge(1) and Si + segmEdge(2)
%lfpEx - external lfp
%rCh - channels to be read
%rawData - boolean variable. If raw data needed rawData = 1(true)
%nlxVer - version ov NLX data (0 - initial, 1 - new, ...)
%
%OUTPUTS
%lfpShft - lfp phased with respect to stimulus moments

rSw = unique(segms(:, 3));%unique sweeps to be read
if rawData %raw data requested
    segmEdge = segmEdge * zavp.rarStep;%left and right shifts from synchro-point (samples)
    segms(:, 1) = segms(:, 1) * zavp.rarStep + 1;%position of synchro-impuls (samples)
else %resampled data requested
    segms(:, 1) = segms(:, 1) + 1;%position of synchro-impuls (samples)
end
segmLen = diff(segmEdge) + 1;%length of needed segments (samples)
segms(:, 1) = round(segms(:, 1));%rounding of noninterg milliseconds

lfpShft = zeros(segmLen, numel(rCh), size(segms, 1));%lfp phased with respect to stimuli moments
for ch = 1:numel(rCh)
    if rawData %need raw data
        lfp = ZavLoadData(zavp.file, hd, rCh(ch), rSw, 0, 'e', nlxVer);%original (raw) data (single channel only)
    else
        lfp = squeeze(lfpEx(:, rCh(ch), rSw));%current channel only
    end
    dataLen = size(lfp, 1);%length of part of data
    for sn = 1:size(segms, 1)
        siPsn = segms(sn, 1);%position of synchro-impuls (ms)
        k = siPsn + segmEdge(1);%negative number equal to number of points to be skipped
        p1 = max(1, k);%first point to be read from lfp
        k = ((1 - k) * (k <= 0)) + 1;%first point in lfpShft for to write
        p2 = min(siPsn + segmEdge(2), dataLen);%last point to be read from lfp
        t = siPsn + segmEdge(2) - dataLen;%negative number equal to number of points to be skipped
        t = t * (t >= 0);%(size(lfpShft, 1) - t) - last point in lfpShft for write
        lfpShft(k:(end - t), ch, sn) = lfp(p1:p2, (rSw == segms(sn, 3)));%lfp phased with respect to stimulus moments
    end
end
