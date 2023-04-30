function [data, hd, ttlIn, spkTS, spkSM] = ZavNrlynx(pf, hd, rCh, nlxVer, strt, stp, irEx)
%[data, ttlIn, hd, spkTS, spkSM] = ZavNrlynx(pf, hd, rCh, nlxVer, strt, stp, irEx)
%read neuralynx. NeuralynxMatlabImportExport_v501 require
%
%INPUTS
%pf - pathname of directory with files *.ncs, *.nev, *.nse (or full pathename of one of file in the directory)
%hd - header (full, for all channels)
%rCh - numbers of channels to be read (one channel only!)
%nlxVer - version ov NLX data (0 - initial, 1 - new, ...)
%strt - start time to read record (ms from beginning of record) - not used
%stp - stop time to read record (ms from beginning of record) - not used
%irEx - exclude interrecords periods if true
%
%OUTPUTS
%data - signal samples (microvoltes)
%hd - file header (information about record)
%ttlIn - moments of synchro-TTL inputs (samples from sweep beginning)
%spkTS - spikes appearence moments (mks)
%spkSM - spikes samples

%define slash
z = find(pf == '\', 1, 'last');%find slash
if isempty(z)
    slashStr = '/';%Linux mode slash
else
    slashStr = '\';%Windows mode slash
end

if (strcmp(pf((end - 3):end), '.ncs') || strcmp(pf((end - 3):end), '.nev') || strcmp(pf((end - 3):end), '.nse'))
    if strcmp(pf((end - 3):end), '.ncs') %ncs file requested
        flNm = pf;%full name specified
    end
    t = find(pf == slashStr, 1, 'last');%last slash
    pf = pf(1:t);%directory paht only
end
if (pf(end) ~= slashStr) %no slash
    pf(end + 1) = slashStr;
end
if ~exist(pf, 'dir') %wrong directory
    %[~, pf] = uigetfile('*.ncs; *.nev; *.nse', 'Select file', pf);%open dialog
    disp(['DIRECTORY ERROR: ', pf])
    return %exit from function ZavNrlynx
end

if (length(rCh) > 1) %more than one channel is requested
    disp('WARNING: only one channel can be read during a single request')
    %за один запрос можно получить данные только одного канала, т.к. в запрос могут попасть каналы с разной частотой дискретизации
    rCh = rCh(1);%only one channel can be read during a single request
end

if isempty(hd) %no header
    rChNum = 0;%counter of channels
    dirCnt = dir(pf);%directory content
    ncsFiles(1:length(dirCnt)) = struct('f', '');%attributes of all channel-files
    lettNumNm = cell(length(dirCnt), 2);%[letter names of channels, numeric name]
    for t = 1:length(dirCnt) %run over files in NLX directory (ncs, we must scan all channels here)
        if ((~dirCnt(t).isdir) && (length(dirCnt(t).name) > 3)) %not directory and name good
            if isequal(dirCnt(t).name(end - 3:end), '.ncs') %data-file of cheetah
                rChNum = rChNum + 1;%increase counter of channels
                recflnmlen = length(dirCnt(t).name);%length of filename
                z = recflnmlen - 4;%character position
                while (z > 1) %run over characters of filename
                    c = double(dirCnt(t).name(z));%double representation of a character
                    if ((c < 48) || (c > 57)) %not a number
                        break;%out of (while z)
                    end
                    z = z - 1;%next character
                end
                ncsFiles(rChNum).f = [pf, dirCnt(t).name];%full name of channel-file
                lettNumNm{rChNum, 1} = dirCnt(t).name(1:z);%letter name of current channel
                lettNumNm{rChNum, 2} = str2double(dirCnt(t).name((z + 1):(recflnmlen - 4)));%numeric name of current channel
                if (isempty(lettNumNm{rChNum, 2}) || ~isfinite(lettNumNm{rChNum, 2}) || any(~isnumeric(lettNumNm{rChNum, 2}))) %not a number
                    disp(['CHANNEL NAME ERROR: ', lettNumNm{rChNum, 2}])
                    return %exit from function ZavNrlynx
                end
            end
        end
    end
    ncsFiles((rChNum + 1):end) = [];%delete empty
    lettNumNm((rChNum + 1):end, :) = [];%delete empty

    %= right order =%
    unqNumNm = unique(vertcat(lettNumNm{:, 2}));%unique numeric names of channels
    if (length(unqNumNm) >= size(lettNumNm, 1)) %no numeric collisions
        [~, rOrd] = sort(vertcat(lettNumNm{:, 2}));%priority of numeric part of names
    else %numeric collisions
        origNumNmOrd = vertcat(lettNumNm{:, 2});%original number names order
        unqLettNm = unique(lettNumNm(:, 1));%unique and sorted letter names of channels
        rOrd = zeros(rChNum, 1);%"right" order of channels
        k = 1;%names counter
        for t = 1:length(unqNumNm) %run over unique number (numeric name is a priority)
            ii = find(origNumNmOrd == unqNumNm(t));%position of the name
            nmPos = length(ii);%number of name position
            if (nmPos == 1) %no collisions
                rOrd(k) = ii;%right order index
                k = k + 1;%next name
            else %number names collision
                jj = zeros(nmPos, 2);%array for order index
                for z = 1:nmPos %run over collided names
                    jj(z, 1) = find(ismember(unqLettNm, lettNumNm(ii(z), 1)));%index of collided names
                end
                [~, jj(:, 2)] = sort(jj(:, 1));%sorted numbers of collided names in array of unique letter names
                rOrd(k + (0:(nmPos - 1))) = ii(jj(:, 2));%right order inex
                k = k + nmPos;%next name
            end
        end
    end
    ncsFiles = ncsFiles(rOrd);%right order of channels
    %= end of (right order) =%
    
    if isempty(rCh) %no channel specified
        rCh = 1;%read the first channels
        for t = 1:length(ncsFiles) %run over files
            if isequal(ncsFiles(t).f, flNm)
                rCh = t;%file detected
                break
            end
        end
    end
    cscHd = Nlx2MatCSC(ncsFiles(rCh).f, [0 0 0 0 0], 1, 1, []);%first and last timestamps and header
    ch_si = (1e6 / GetNlxHeadParametr(cscHd, 'SamplingFrequency'));%shortest sampling interval (mks)
    adBitVolts = GetNlxHeadParametr(cscHd, 'ADBitVolts');%multiplier to convert from samples to volts (lfp)
    dspDelay_mks = GetNlxHeadParametr(cscHd, 'DspFilterDelay_µs');%DspFilterDelay_µs (lfp)
    if ~isempty(dspDelay_mks)
        dspDelay_mks = dspDelay_mks * double(isequal(GetNlxHeadParametr(cscHd, 'DspDelayCompensation'), 'Disabled'));%DspDelayCompensation (lfp)
    end
    inpInvert = double(strcmp(GetNlxHeadParametr(cscHd, 'InputInverted'), 'True'));%input inverted
else %we have a header
    ncsFiles(1:hd.nADCNumChannels) = struct('f', '');%attributes of all channel-files
    for ch = 1:hd.nADCNumChannels %run over channels
        ncsFiles(ch).f = [pf, hd.recChNames{ch}, '.ncs'];%full names of channel-files
    end
    if isempty(rCh) %no channel specified
        rCh = 1;%read the first channels
        for t = 1:length(ncsFiles) %run over files
            if isequal(ncsFiles(t).f, flNm)
                rCh = t;%file detected
                break
            end
        end
    end
    ch_si = hd.ch_si(rCh);%sampling interval of the requested channel (mks)
    adBitVolts = hd.adBitVolts(rCh);%multiplier to convert from samples to volts (lfp)
    dspDelay_mks = hd.dspDelay_mks(rCh);%DspFilterDelay_µs (lfp)
    inpInvert = hd.inverted(rCh);%input inverted
end
rChNum = length(rCh);%number of requested channels

fileToRead = [pf, 'Events.nev'];%pathname of file to be read
if exist(fileToRead, 'file') %events-file exist
    [evntTStmp, ttl, evntStr] = Nlx2MatEV(fileToRead, [1 0 1 0 1], 0, 1, []);%read events, timestamps and event strings
end
    
%= read data from ncs-files =%
cscTStmp = Nlx2MatCSC(ncsFiles(rCh).f, [1 0 0 0 0], 0, 1, []);%timestamps of recorded samples (read largest file)
if any(diff(cscTStmp, 2) > 5) %datasections lost
    disp('large error of timestamp')%datasections lost
end
if (isempty(strt) || (strt < 0)) %default value
    strt = 0;%begin read from recordation start (ms from record begin)
end
if (isempty(stp) || isequal(stp, 'e') || (stp > ((diff(cscTStmp([1, end])) + (512 * ch_si)) / 1e3))) %read upto end of recordation
    stp = (diff(cscTStmp([1, end])) + (512 * ch_si)) / 1e3;%read until end of recordation (ms)
end
rRecrd = [0, 0];%number of 512-samples-records to be read (initiation of array)

fileToRead = ncsFiles(rCh).f;%[pf, 'CSC', num2str(ch), '.ncs'];%pathname of file to be read
if exist(fileToRead, 'file') %requested ncs-file exist
    if (cscTStmp(1) > 0) %we can read file fully
        %= get requested samples only =%
        rRecrd(1) = find(cscTStmp <= (cscTStmp(1) + (strt * 1e3)), 1, 'last');%number of first 512-samples-record to be read
        rRecrd(2) = find(cscTStmp <= (cscTStmp(1) + (stp * 1e3)), 1, 'last');%number of last 512-samples-record to be read
        [nValSampl, smpl] = Nlx2MatCSC(fileToRead, [0 0 0 1 1], 0, 2, rRecrd);%get requested samples
        %= end of get (requested samples only) =%
        % %smpl = Nlx2MatCSC(fileToRead, [0 0 0 0 1], 0, 1, []);%get all samples

        if any(diff(cscTStmp, 2) > 5) %%(2)complex way (segment-mode or datasections lost; use timestamps)
            %disp('large error of timestamp')
            mStStRec = FindStStStemp(evntStr, evntTStmp, cscTStmp, ch_si);%find moments of "Starting Recording" and "Stopping Recording" (microseconds)
            copy_cscTStmp = cscTStmp(rRecrd(1):rRecrd(2));%copy
            %= exclude interrecords periods =%
            if irEx %exclude interrecords periods if true
                for z = (numel(mStStRec) - 1):-2:2 %run over start-stop events
                    jj = (copy_cscTStmp >= mStStRec(z));%number of timestamps satisfying conditions
                    copy_cscTStmp(jj) = copy_cscTStmp(jj) - (mStStRec(z) - mStStRec(z - 1)) + (512 * ch_si);%exclude interrecords periods
                end
            end
            %=end of (exclude interrecords periods) =%

            copy_cscTStmp = copy_cscTStmp - copy_cscTStmp(1);%adduction the first sample to zeros

            data(1:floor(diff(copy_cscTStmp([1, end])) / ch_si), 1) = NaN;%memory preallocation
            for t = 1:length(copy_cscTStmp) %run over timestamps
                k = floor(copy_cscTStmp(t) / ch_si) + 1;
                data(k + (0:(nValSampl(t) - 1)), 1) = smpl(1:nValSampl(t), t);%samples
            end
            data((k + 512):end, 1) = [];%delete excess
        else %(1)simple way (concatenation)
            data = smpl(:);%samples
        end

        data = data * adBitVolts * 1e6;%microvoltes (lfp)
        if (inpInvert >= 1) %inverted signal
            data = -1 * data;%back inverse
        end
        
        %cut nose and tail from data
        cutSmplN = round(((cscTStmp(1) + (strt * 1e3)) - cscTStmp(rRecrd(1))) / ch_si) - 1;%number of nose-samples to be cutted
        if (cutSmplN > 0) %need to cut nose
            data(1:cutSmplN, :) = [];%cut nose
        end
        cutSmplN = 512 - round(((cscTStmp(1) + (stp * 1e3)) - cscTStmp(rRecrd(2))) / ch_si) - 1;%number of tail-samples to be cutted
        if (cutSmplN > 0) %need to cut tail
            data((end - cutSmplN):end, :) = [];%cut tail
        end
    else %empty record
        data = [];%no data exist
    end
else %file not found
    data = [];%no data exist
end
%= end of (read data from ncs-files) =%

if ((nargout > 1) && isempty(hd)) %header requested
    hd.fFileSignature = 'Neuralynx';
    hd.nOperationMode = 3;%data were acquired in gap-free mode (continuous record)
    hd.lActualEpisodes = 1;%number of sweeps (for compatibility with abfload)
    hd.nADCNumChannels = length(ncsFiles);%total number of channels
    hd.adBitVolts = zeros(hd.nADCNumChannels, 1);%multiplier to convert from samples to volts (lfp)
    hd.dspDelay_mks = zeros(hd.nADCNumChannels, 1);%DspFilterDelay_µs (lfp)
    hd.adBitVoltsSpk = zeros(hd.nADCNumChannels, 1);%multiplier to convert from samples to volts (spikes)
    hd.dspDelay_mksSpk = zeros(hd.nADCNumChannels, 1);%DspFilterDelay_µs (spikes)
    hd.alignmentPt = zeros(hd.nADCNumChannels, 1);%spike samples back (from peak, including peak point)
    hd.inverted = zeros(hd.nADCNumChannels, 1);%input inverted
    hd.recChUnits = cell(hd.nADCNumChannels, 1);%mesurement units
    hd.recChNames = cell(hd.nADCNumChannels, 1);%name of channels
    hd.recChDurations = zeros(hd.nADCNumChannels, 1);%duration of the record (microseconds)
    hd.ch_si = NaN(hd.nADCNumChannels, 1);%sample interval (mks)
    hd.dataPtsPerChan = zeros(hd.nADCNumChannels, 1);%samples per channel
    hd.DSPLowCutFilterEnabled = true(hd.nADCNumChannels, 1);%low-cut filter enabled
    hd.DspLowCutFrequency = zeros(hd.nADCNumChannels, 1);%low-cut filter frequency
    hd.DSPHighCutFilterEnabled = true(hd.nADCNumChannels, 1);%high-cut filter enabled
    hd.DspHighCutFrequency = zeros(hd.nADCNumChannels, 1);%high-cut filter frequency

    %read headers
    for ch = 1:hd.nADCNumChannels %run over channels
        %%% CSC(ncs)-files (lfp) %%%
        fileToRead = ncsFiles(ch).f;%[pf, 'CSC', num2str(ch), '.ncs'];%pathname of file to be read
        if exist(fileToRead, 'file') %requested file with lfp exist
            [cscTStmp, cscHd] = Nlx2MatCSC(fileToRead, [1 0 0 0 0], 1, 1, []);%first and last timestamps and header
            if (cscTStmp(1) > 0) %existing recoredation
                hd.ch_si(ch) = (1e6 / GetNlxHeadParametr(cscHd, 'SamplingFrequency'));%sample interval (mks)
                hd.dataPtsPerChan(ch) = length(cscTStmp) * 512;%samples per channel
                hd.recChDurations(ch) = hd.dataPtsPerChan(ch) * hd.ch_si(ch);%diff(cscTStmp([1, end]));%duration of the record (microseconds)
            end
            hd.adBitVolts(ch) = GetNlxHeadParametr(cscHd, 'ADBitVolts');%multiplier to convert from samples to volts (lfp)
            dspDelay_mks = GetNlxHeadParametr(cscHd, 'DspFilterDelay_µs');%DspFilterDelay_µs (lfp)
            if isempty(dspDelay_mks)
                dspDelay_mks = Inf;
            end
            hd.dspDelay_mks(ch) = dspDelay_mks;%DspFilterDelay_µs (lfp)
            hd.dspDelay_mks(ch) = hd.dspDelay_mks(ch) * double(isequal(GetNlxHeadParametr(cscHd, 'DspDelayCompensation'), 'Disabled'));%DspDelayCompensation (lfp)
            hd.recChUnits{ch} = 'µV';%mesurement units
            hd.recChNames{ch} = GetNlxHeadParametr(cscHd, 'AcqEntName');%name of channels
            z = find(fileToRead == '\');
            if ~isequal(fileToRead((z(end) + 1):(end - 4)), hd.recChNames{ch}) %manually renamed file
                hd.recChNames{ch} = fileToRead((z(end) + 1):(end - 4));%for manually renamed files
            end
            hd.inverted(ch) = double(strcmp(GetNlxHeadParametr(cscHd, 'InputInverted'), 'True'));%input inverted
            hd.DSPLowCutFilterEnabled(ch) = isequal(GetNlxHeadParametr(cscHd, 'DSPLowCutFilterEnabled'), 'True');%low-cut filter enabled
            hd.DspLowCutFrequency(ch) = GetNlxHeadParametr(cscHd, 'DspLowCutFrequency');%low-cut filter frequency
            hd.DSPHighCutFilterEnabled(ch) = isequal(GetNlxHeadParametr(cscHd, 'DSPHighCutFilterEnabled'), 'True');%low-cut filter enabled
            hd.DspHighCutFrequency(ch) = GetNlxHeadParametr(cscHd, 'DspHighCutFrequency');%low-cut filter frequency
        end

        %%% SE(nse)-files (spikes) %%%
        fileToRead = [pf, 'SE', num2str(ch), '.nse'];%pathname of file to be read
        if exist(fileToRead, 'file') %requested file with spikes exist
            spkHd = Nlx2MatSpike(fileToRead, [0 0 0 0 0], 1, 1, []);%header (spikes)
            hd.adBitVoltsSpk(ch) = GetNlxHeadParametr(spkHd, 'ADBitVolts');%multiplier to convert from samples to volts (spikes)
            hd.dspDelay_mksSpk(ch) = GetNlxHeadParametr(spkHd, 'DspFilterDelay_µs');%DspFilterDelay_µs (spikes)
            hd.dspDelay_mksSpk(ch) = hd.dspDelay_mksSpk(ch) * double(isequal(GetNlxHeadParametr(spkHd, 'DspDelayCompensation'), 'Disabled'));%DspDelayCompensation (spikes)
            hd.alignmentPt(ch) = GetNlxHeadParametr(spkHd, 'AlignmentPt');%spike samples back (from peak, including peak point)
        end
    end
    hd.dataPts = sum(hd.dataPtsPerChan);%total number of recorded samples
    hd.si = min(hd.ch_si(hd.recChDurations > 0));%shortest sampling interval (mks)
    hd.fADCSampleInterval = hd.si;%sample interval (mks)

    if exist([pf, 'Events.nev'], 'file') %events-file exist
        hd.TTLs = ttl;%ttl events
        hd.EventStrings = evntStr;%text of events
    end
    
    %log-file read
    fid = fopen([pf, 'CheetahLogFile.txt']);%open cheetah log-file
    logy = fread(fid, 'char');%read log-file
    logy = char(logy');
    fclose(fid);

    jj = strfind(logy, '*-');
    tmS = cell(numel(jj), 3);
    if (nlxVer == 0) %initial version of cheetah logs
        for t = 1:numel(jj)
            z = strfind(logy(jj(t):(jj(t) + 50)), ' - ');
            bufStr = textscan(logy((jj(t) + 2):(jj(t) + z(2) - 2)), '%s%s%s');
            tmS(t, 1) = bufStr{1}(1);%time of current event (hh:mm:ss.ms)
            tmS{t, 2} = str2double(bufStr{3}{1});%timestamp of current event (microseconds from Neuralynx on)
            tmS{t, 3} = logy((jj(t) + z(2) + 2):(jj(t) + z(2) + 2 + 33));%current event string
        end
        ymdhms = datevec(tmS{1, 1}, 'HH:MM:SS.FFF');%date as numels
    else %next version of cheetah logs
        for t = 1:numel(jj)
            z = strfind(logy(jj(t):(jj(t) + 50)), ' - ');
            tmS{t, 1} = logy((jj(t) + 2):(jj(t) + z(1) - 2));%date and time of current event (yyyy/mm/dd HH:MM:SS)
            tmS{t, 2} = str2double(logy((jj(t) + 2 + z(1)):(jj(t) + z(2) - 2)));%timestamp of current event (microseconds from Neuralynx on)
            tmS{t, 3} = logy((jj(t) + z(2) + 2):(jj(t) + z(2) + 2 + 33));%current event string
        end
        ymdhms = datevec(tmS{1, 1}, 'yyyy/mm/dd HH:MM:SS');%date as numels
    end
    nlxBeginS = ymdhms(4) * 3600 + ymdhms(5) * 60 + ymdhms(6);%time of Neuralynx turned ON (seconds from day begin)
    
    n = 0;%number of first event 'Starting Recording'
    for t = 1:length(evntStr)
        if isequal(evntStr{t}, 'Starting Recording')
            n = t;%number of first event 'Starting Recording'
            break;%out of "for t"
        end
    end
    [~, t] = min(abs(vertcat(tmS{:, 2}) - evntTStmp(n)));%nearest timestamp
    jj = false;
    z = -1;
    while (~jj)
        z = z + 1;
        jj = strcmp(tmS{t + z, 3}, 'AcquisitionControl::StartRecording');%go forward to find
        if ~jj
            jj = strcmp(tmS{t - z, 3}, 'AcquisitionControl::StartRecording');%go backward to find
            if jj
                z = -z;
            end
        end
    end
    t = t + z;%number of row with pure time of start
    if (nlxVer == 0) %initial version of cheetah logs
        rcTm = datenum(['2011:01:01 ', tmS{t, 1}], 'yyyy:mm:dd HH:MM:SS.FFF') - datenum('2011:01:01 00:00:00.000', 'yyyy:mm:dd HH:MM:SS.FFF');%days after day begin
        hd.recTime(1) = rcTm * (24 * 60 * 60);%record stop time (seconds after day begin)
        if ((hd.recTime(1) - nlxBeginS) <= 0) %new day began
            hd.recTime(1) = hd.recTime(1) + (24 * 60 * 60);%record stop time (seconds after day begin)
        end
        %hd.recTime(1) = (tmS{t, 2}  / 1e6) + nlxBeginMKS;%record start time (seconds from day begin (seconds after day begin))
    else
        rcTm = datenum(tmS{t, 1}, 'yyyy/mm/dd HH:MM:SS') - datenum([tmS{1, 1}(1:11), ' 00:00:00'], 'yyyy/mm/dd HH:MM:SS');%days after day begin
        hd.recTime(1) = rcTm * (24 * 60 * 60);%record stop time (seconds after day begin)
    end
    
    n = 0;%number of last event 'Stopping Recording'
    for t = length(evntStr):-1:1
        if isequal(evntStr{t}, 'Stopping Recording')
            n = t;%number of last event 'Stopping Recording'
            break;%out of "for t"
        end
    end
    [~, t] = min(abs(vertcat(tmS{:, 2}) - evntTStmp(n)));%nearest timestamp
    jj = false;
    z = -1;
    while (~jj)
        z = z + 1;
        jj = strcmp(tmS{t + z, 3}, 'AcquisitionControl::StopRecording(');
        if ~jj
            jj = strcmp(tmS{t - z, 3}, 'AcquisitionControl::StopRecording(');
            if jj
                z = -z;
            end
        end
    end
    t = t + z;%number of row with pure time of stop
    if (nlxVer == 0) %initial version of cheetah logs
        rcTm = datenum(['2011:01:01 ', tmS{t, 1}], 'yyyy:mm:dd HH:MM:SS.FFF') - datenum('2011:01:01 00:00:00.000', 'yyyy:mm:dd HH:MM:SS.FFF');%days after day begin
        hd.recTime(2) = rcTm * (24 * 60 * 60);%record stop time (seconds after day begin)
        if ((hd.recTime(2) - nlxBeginS) <= 0) %new day began
            hd.recTime(2) = hd.recTime(2) + (24 * 60 * 60);%record stop time (seconds after day begin)
        end
    else
        rcTm = datenum(tmS{t, 1}, 'yyyy/mm/dd HH:MM:SS') - datenum([tmS{1, 1}(1:11), ' 00:00:00'], 'yyyy/mm/dd HH:MM:SS');%days after day begin
        hd.recTime(2) = rcTm * (24 * 60 * 60);%record stop time (seconds after day begin) (not good precision)
    end
%     hd.recTime(2) = hd.recTime(1) + diff(mStStRec([1, end])) / 1e6;%time of record stop (better precision)
%     [~, z] = max(hd.hd.dataPtsPerChan);%longest recordation
%     hd.recTime(2) = hd.recTime(1) + hd.dataPtsPerChan(z) * hd.ch_si(z) / 1e6;%time of record stop (the best precision)
end %end of (if ((nargout > 1) && isempty(hd)) %header requested)


if (nargout > 2) %synchro-signals requested
    frstTStmp = Inf;%first timestamp of the record
    for ch = 1:length(ncsFiles) %run over ncs-files
        if (hd.recChDurations(ch) > 0) %record exist
            cscTStmp = Nlx2MatCSC(ncsFiles(ch).f, [1 0 0 0 0], 0, 3, 1);%first and last timestamps and header
            if ((cscTStmp(1) > 0) && (cscTStmp(1) < frstTStmp)) %valid timestamps
                frstTStmp = cscTStmp(1);%first timestamp of the record
            end
        end
    end
    ch = find(hd.ch_si == hd.si, 1, 'first');%channel with minimal sampling interval
    cscTStmp = Nlx2MatCSC(ncsFiles(ch).f, [1 0 0 0 0], 0, 1, []);%timestamps of recorded samples (read largest file)
    mStStRec = FindStStStemp(evntStr, evntTStmp, cscTStmp, hd.si);%find moments of "Starting Recording" and "Stopping Recording" (microseconds)
    if (frstTStmp > 0) %good record
        inEvntOn = zeros(1, length(evntStr));%numbers of events when input ports changed
        for t = 1:length(evntStr) %run over event strings
            inEvntOn(t) = t * double(~isempty(strfind(evntStr{t}, 'Input')));%find input events
            %inEvntOn(t) = t * double(~isempty(strfind(evntStr{t}, 'Output')));%find output events
        end
        inEvntOff = inEvntOn((ttl <= 0) & (inEvntOn > 0));%number of events when input TTL ports are in state 'OFF'
        inEvntOn = inEvntOn((ttl > 0) & (inEvntOn > 0));%number of events when input TTL ports are in state 'ON'

        ttlPrtOn = unique(evntStr(inEvntOn));%different TTL ports (input ports only) in state 'ON'
        ttlIn(1:length(ttlPrtOn)) = struct('t', zeros(numel(inEvntOn), 2));%initialization
        origTTL(1:length(ttlPrtOn)) = struct('t', []);%initialization
        for ch = 1:length(ttlPrtOn) %run over different TTL ports
            z = 1;%counter of synchroimpulses
            for t = inEvntOn %run over inputs events when TTL set to On'
                if strcmp(evntStr{t}, ttlPrtOn{ch}) %right number of inputs port
                    ttlIn(ch).t(z, 1) = evntTStmp(t);%"on" (allStims(:, 1)) stimulus (mks from NLX switch on)
                    for n = inEvntOff(inEvntOff > t) %run over input events when TTL set to 'Off'
                        if strcmp(evntStr{n}(1:(end - 12)), ttlPrtOn{ch}(1:(end - 12))) %right number of inputs port
                            ttlIn(ch).t(z, 2) = evntTStmp(n);%"off"(allStims(:, 2)) stimulus (mks from NLX switch on)
                            z = z + 1;%counter of synchroimpulses
                            break;%out of (for n)
                        end
                    end
                end
            end
            ttlIn(ch).t(z:end, :) = [];%delete excess
            origTTL(ch).t = ttlIn(ch).t;%original timestamps of input TTLs without interrecords periods exclusion

            %= exclude interrecords periods =%
            if irEx %exclude interrecords periods if true
                for z = (numel(mStStRec) - 1):-2:2 %run over start-stop events
                    jj = (ttlIn(ch).t(:, 1) >= mStStRec(z));%number of timestamps satisfying conditions
                    ttlIn(ch).t(jj, :) = ttlIn(ch).t(jj, :) - (mStStRec(z) - mStStRec(z - 1)) + (512 * hd.si);%stimulus moments from record begin
                end
            end
            %= end of (exclude interrecords periods) =%

            ttlIn(ch).t = ttlIn(ch).t - frstTStmp;%adduction to zeros (first sample, mks from start-record-event)
            ttlIn(ch).t = ttlIn(ch).t / hd.si;%convert stimulus moments to samples from record begin (from 0 because of arbitrary position of synchro-event in relation to sampling time)
            origTTL(ch).t = origTTL(ch).t - frstTStmp;%adduction to zeros (first sample, mks from start-record-event)
        end
        hd.inTTL_timestamps = origTTL;%original timestamps of input TTLs (adducted to zeros - first sample is zero)
        hd.sweepStartInPts = (mStStRec - mStStRec(1)) / hd.si;%the start times of sweeps in sample points (from beginning of recording)
    else
        ttlIn = [];%no events-file exist
    end
end %end of (if (nargout > 2) %synchro-signals requested)


if (nargout > 3) %timestamps of spikes is requested
    spkTS(1:rChNum, 1) = struct('tStamp', []);%spikes moment (mks)
    fileToRead = [pf, 'SE', num2str(ch), '.nse'];%pathname of file to be read
    if exist(fileToRead, 'file') %requested file with spikes exist
        if isempty(hd) %no header
            spkHd = Nlx2MatSpike(fileToRead, [0 0 0 0 0], 1, 1, []);%header (spikes)
            adBitVoltsSpk = GetNlxHeadParametr(spkHd, 'ADBitVolts');%multiplier to convert from samples to volts (spikes)
            dspDelay_mksSpk = GetNlxHeadParametr(spkHd, 'DspFilterDelay_µs');%DspFilterDelay_µs (spikes)
            dspDelay_mksSpk = dspDelay_mksSpk * double(isequal(GetNlxHeadParametr(spkHd, 'DspDelayCompensation'), 'Disabled'));%DspDelayCompensation (spikes)
        else %we have a header
            adBitVoltsSpk = hd.adBitVoltsSpk(ch);%multiplier to convert from samples to volts (spikes)
            dspDelay_mksSpk = hd.dspDelay_mksSpk(ch);%DspFilterDelay_µs (spikes)
            dspDelay_mks = hd.dspDelay_mks(ch);%DspDelayCompensation (scs)
        end
        spkTmStmp = Nlx2MatSpike(fileToRead, [1 0 0 0 0], 0, 1, []);%timestamps of spikes
%             spkTmStmp = spkTmStmp((spkTmStmp >= rStmps(1)) & (spkTmStmp <= rStmps(2)));%wanted spikes only
%                 %or:
%                 %spkTmStmp = Nlx2MatSpike(fileToRead, [1 0 0 0 0], 0, 4, rStmps);

        spkTS(n).tStamp = spkTmStmp - evntTStmp(1) - round((dspDelay_mksSpk + dspDelay_mks) / 2);%mks from record start
    end
    if (nargout > 4) %spikes time course requested
        spkSM(1:rChNum, 1) = struct('shape', []);%spikes samples
%         spkSM(n).shape = squeeze(Nlx2MatSpike(fileToRead, [0 0 0 0 1], 0, 4, rStmps));
        spkSM(n).shape = squeeze(Nlx2MatSpike(fileToRead, [0 0 0 0 1], 0, 1, []));
        spkSM(n).shape = spkSM(n).shape * adBitVoltsSpk * 1e6;%samples to microvolts
    end
end


function mStStRec = FindStStStemp(evntStr, ststTStmp, cscTS, si)
%find moments of "Starting Recording" and "Stopping Recording"
%
%INPUTS
%evntStr - strings with events description
%ststTStmp - timestampes of start-stop and ttl events
%cscTS - timestamps of recorded samples (mks, largest file)
%si - sample interval (mks)
%
%OUTPUTS
%mStStRec - numbers of timestamp of start and stop recordings
%

mStStRec = zeros(length(evntStr), 1);%memory preallocation
z = 1;
for t = 1:length(evntStr)
    if strcmp(evntStr{t}, 'Starting Recording')
        mStStRec(z) = ststTStmp(t);%start timestamp (number)
        z = z + 1;
    end
    if strcmp(evntStr{t}, 'Stopping Recording')
        mStStRec(z) = ststTStmp(t);%stop timestamp (number)
        z = z + 1;
    end
end
mStStRec(z:end, :) = [];%delete excess
for t = 1:2:length(mStStRec) %run over start timestamps
    [~, z] = min(abs(cscTS - mStStRec(t)));%number of nearest timestemp of samples-block
    mStStRec(t) = cscTS(z);%nearest timestemp of samples-block
end
for t = 2:2:length(mStStRec) %run over stop timestamps
    [~, z] = min(abs(cscTS - mStStRec(t)));%number of nearest timestemp of samples-block
    mStStRec(t) = cscTS(z) + 512 * si;%nearest timestemp of samples-block
end

