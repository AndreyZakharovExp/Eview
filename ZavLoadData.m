function [data, hd, synchro, spkTS, spkSM] = ZavLoadData(flNm, inHd, rCh, rSw, strt, stp, nlxVer)
%[data, hd, synchro, spkTS, spkSM] = ZavLoadData(flNm, inHd, rCh, rSw, strt, stp, nlxVer)
%unified interface for load files (abf-files, Neuralynx, edf+)
%
%INPUTS
%flNm - pathname of file
%inHd - input header
%rCh - numbers of channels to be read
%rSw - numbers of sweeps to be read
%strt - start time (ms from record begin)
%stp - stop time (ms from record begin)
%nlxVer - version ov NLX data (0 - initial, 1 - new, ...)
%
%OUTPUTS
%data - signal samples
%hd - file header (information about record)
%synchro - stimulus moments (samples from sweep beginning)
%spkTS - spikes appearence moments (mks)
%spkSM - spikes samples (Neuralynx)
%

t = find(flNm(1:end) == '.', 1, 'last');
if (~isempty(t) && ((length(flNm) - t) < 7))
    ext = flNm(t:end);%file extention exist
else
    ext = '';%no file extention
end

if (nargout > 1) %synchro requested
    synchro = struct('t', []);%initial empty
end
if (nargout > 3) %spikes timestampes requested
    spkTS = [];%initial empty
end
if (nargout > 4) %spikes samples requested
    spkSM = [];%initial empty
end

if isempty(strt) %no start time specified
    strt = 0;%begin read from recordation start
end
if isempty(stp) %no stop time specified
    stp = 'e';%read until end of recordation
end
if (strt < 0) %before record
    strt = 0;%first sample time (ms)
end
if ~isempty(inHd) %header exist
    if ischar(rCh)
        ch = 1;%read first channel
    else
        ch = rCh(1);%read first channel
    end
    if (stp > (inHd.recChDurations(ch) * 1e-3)) %after record
        stp = inHd.recChDurations(ch) * 1e-3;%last sample time (ms)
    end
    if (strt > (inHd.recChDurations(ch) * 1e-3)) %wrong start-stop time
        disp('start-stop time error')
        return %out of ZavLoadData
    end
end
if ((stp < 0) || (strt > stp)) %wrong start-stop time
    disp('start-stop time error')
    return %out of ZavLoadData
end

switch ext
    case '.abf' %axon binary file
        if isempty(rCh) %no channel number specified
            rCh = 'a';%read all channels
        end
        if isempty(rSw) %no sweep number specified
            rSw = 'a';%read all sweeps
        end
        strt = strt / 1e3;%convert to seconds
        if (stp ~= 'e')
            stp = stp / 1e3;%convert to seconds
        end
        if isempty(inHd) %no header specified
            [~, ~, hd] = abfload(flNm, 'channels', 'a', 'sweeps', 1, 'start', 0, 'stop', 'e');%get header
        else
            hd = inHd;%same header
        end
        if (~isequal(rCh, 'a') && isnumeric(rCh))
            abfCh = cell(1, length(rCh));%names of requested channels
            for t = 1:length(rCh)
                abfCh{t} = hd.recChNames{rCh(t)};%names of requested channels
            end
        else %rCh containes channel names (must be cell or string) or string 'a'
            abfCh = rCh;%names of rquested channels
        end
        if (hd.nOperationMode == 3) %gap-free mode
            %what if ((strt ~= 0) || (stp ~= 'e')) ?
            data = abfload(flNm, 'channels', abfCh, 'sweeps', rSw, 'start', strt, 'stop', stp);
            hd.lActualEpisodes = 1;%number of sweeps
        elseif any(hd.nOperationMode == [1, 2, 5]) %event-driven mode
            data = abfload(flNm, 'channels', abfCh, 'sweeps', rSw, 'start', 0, 'stop', 'e');
            recLen = size(data, 1);%total samples number
            pL = round(strt * 1e6 / hd.si) + 1;%first sample
            if (stp == 'e') %read up to end
                pR = recLen;%final sample
            else %read partially
                pR = round(stp * 1e6 / hd.si);%final sample
            end
            data = data(pL:pR, :, :);%cut requested time
        end
        hd.ch_si(1:hd.nADCNumChannels) = hd.si;%same sampling rate on each channels
        hd.recChDurations = (hd.dataPtsPerChan .* ones(hd.nADCNumChannels, 1)) * hd.si;%duration of the record (microseconds)
        hd.dataPtsPerChan = hd.dataPtsPerChan .* ones(hd.nADCNumChannels, 1);%samples per channel
    case {'.ncs', '.nse', '.nev', ''} %Neuralynx
        if ischar(rCh) %string value (for example 'a' like for abf)
            disp('rCh is not valid')
            return %out of ZavLoadData
        end
        if isempty(nlxVer) %version of Neuralynx data
            nlxVer = 0;%old version of Neuralynx Cheetah
        end
        switch nargout
            case 5 %all variables requested
                [data, hd, synchro, spkTS, spkSM] = ZavNrlynx(flNm, inHd, rCh(1), nlxVer, strt, stp, false);%load 
            case 4 %first 4 variables requested
                [data, hd, synchro, spkTS] = ZavNrlynx(flNm, inHd, rCh(1), nlxVer, strt, stp, false);%load 
            case 3 %first 3 variables requested
                [data, hd, synchro] = ZavNrlynx(flNm, inHd, rCh(1), nlxVer, strt, stp, false);%load 
            case {2, 1} %first 2 variables or samples only requested
                [data, hd] = ZavNrlynx(flNm, inHd, rCh(1), nlxVer, strt, stp, false);%load 
        end
        if (numel(rCh) > 1) %WARNING: only one channel can be read during a single request
            data = [data, zeros(size(data, 1), numel(rCh) - 1)];%only one channel can be read during a single request
            for ch = 2:length(rCh) %run over channels requested
                smpl = ZavNrlynx(flNm, hd, rCh(ch), nlxVer, strt, stp, false);%load 
                %разные каналы могут содержать разное количество отсчётов (например, из-за разной частоты дискретизации)
                data(:, ch) = smpl;%only one channel can be read during a single request
            end
        end
    case '.daq' %daq-files;
        switch nargout
            case {3, 4, 5}
                [data, hd] = ZavDaqRead(flNm, inHd, rCh, rSw, strt, stp);%read specified channels
            case {1, 2}
                data = ZavDaqRead(flNm, inHd, rCh, rSw, strt, stp);%read specified channels
        end
    case '.edf' %eeg in "european data format"
        [data, hd] = ZavManEEGload(flNm, rCh);%load EEG in "european data format"
    case '.oebin' %OpenEphys files
        [data, hd] = ZavReadOEphys(flNm, inHd, rCh, [], [], []);%load OpenEphys data
    otherwise
        disp('file type error')
end
