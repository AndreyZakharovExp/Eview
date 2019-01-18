function varargout = Eview(varargin)
% EVIEW MATLAB code for Eview.fig
%      EVIEW, by itself, creates a new EVIEW or raises the existing
%      singleton*.
%
%      H = EVIEW returns the handle to a new EVIEW or the handle to
%      the existing singleton*.
%
%      EVIEW('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EVIEW.M with the given input arguments.
%
%      EVIEW('Property','Value',...) creates a new EVIEW or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Eview_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Eview_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Eview

% Last Modified by GUIDE v2.5 17-Jan-2019 18:54:47

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Eview_OpeningFcn, ...
                   'gui_OutputFcn',  @Eview_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before Eview is made visible.
function Eview_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Eview (see VARARGIN)

% Choose default command line output for Eview
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Eview wait for user response (see UIRESUME)
% uiwait(handles.Eview);


% --- Outputs from this function are returned to the command line.
function varargout = Eview_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
evi = handles.Eview;%main form handle
handles.LFPscb = annotation('textarrow', 'Position', [0.025, 0.15 0, 0.1], 'LineWidth', 5, 'Units', 'normalized', 'TextRotation', 90, 'HeadStyle', 'none',...
    'String', '0 \muV', 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');%LFP-scale bar
handles.MUAscb = annotation('textarrow', 'Position', [0.55, 0.15 0, 0.1], 'LineWidth', 5, 'Units', 'normalized', 'TextRotation', 90, 'HeadStyle', 'none',...
    'String', '0 units', 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');%MUA-scale bar
setappdata(evi, 'dataCurs', []);

handles.CSD_clb = colorbar('peer', handles.LFP_ax, 'Location', 'manual', 'Units', 'normalized', ... %CSD-colorbar
    'Position', [0.027, 0.686, 0.01, 0.17], 'YAxisLocation', 'left', 'FontSize', 6);%CSD-colorbar
handles.Spk_clb = colorbar('peer', handles.Spk_ax, 'Location', 'manual', 'Units', 'normalized', ... %MUA-colorbar
    'Position', [0.527, 0.686, 0.01, 0.17], 'YAxisLocation', 'left', 'FontSize', 6);%MUA-colorbar

set(evi, 'WindowButtonDownFcn', {@Eview_WindowButtonDownFcn, handles})%reaction on button down
set(evi, 'WindowButtonMotionFcn', {@Eview_WindowButtonMotionFcn, handles})%reaction on button motion
set(evi, 'WindowButtonUpFcn', {@Eview_WindowButtonUpFcn, handles})%reaction on button up

guidata(hObject, handles);%update handles structure


% --------------------------------------------------------------------
function OpenData_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to OpenData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

evi = handles.Eview;%main form handle
%nts = handles.NotesForUser;%notes text
pth = getappdata(evi, 'pth');%directory
flNm = getappdata(evi, 'flNm');%filename

if isappdata(evi, 'brst')
    ynBtn = questdlg('Save resultes?', 'Save on exit', 'Yes', 'No', 'No');%ask what to do
    if isequal(ynBtn, 'Yes') %save requested
        SaveBursts_ClickedCallback([], [], handles);%save bursts
    end
end
        
[flNm, pth] = uigetfile({'*.mat', 'mat'}, 'Select file', [pth, flNm, '.mat']);%open dialog
if (isequal(flNm, 0) && isequal(pth, 0)) %no file choosed
    %set(nts, 'String', 'no file choosed')
    return
end

z = find(flNm == '.', 1, 'last');
if ~isequal(flNm(z:end), '.mat')
    %set(nts, 'String', 'unknown file type')
    return
end
flNm = flNm(1:(z - 1));%delete expansion

% ClearData(handles);%clear application content
setappdata(evi, 'pth', pth);%directory
setappdata(evi, 'flNm', flNm);%filename

varInMatFile = who(matfile([pth, flNm, '.mat']));
if ismember('zp', varInMatFile)
    load([pth, flNm, '.mat'], 'lfp', 'spks', 'hd', 'zp', 'lfpVar')
    zavp = zp;
else
    load([pth, flNm, '.mat'], 'lfp', 'spks', 'hd', 'zavp', 'lfpVar')
end
set(evi, 'Name', ['Eview - ', flNm])
        
%lfp = double(lfp);
prg = str2num(get(handles.SpkPrg, 'String'));
spksOrig = spks;
for sw = 1:hd.lActualEpisodes
    for ch = 1:hd.nADCNumChannels
        if isfield(spksOrig, 's')
            spksOrig(ch, sw).tStamp = spksOrig(ch, sw).s;%
            spks(ch, sw).tStamp = spks(ch, sw).s;%
        end
        if isfield(spksOrig, 'ampl')
            ii = double(spksOrig(ch, sw).ampl) <= (-lfpVar(ch) * prg);
            spks(ch, sw).tStamp = spksOrig(ch, sw).tStamp(ii);
            spks(ch, sw).ampl = spksOrig(ch, sw).ampl(ii);
        end
    end
end

%segms(:, 1) - position of synchro-event (ms)
%segms(:, 2) - number of channel where syncro-event was detected
%segms(:, 3) - number of sweeps where syncro-event was detected
%segms(:, 4) - number of stimulus (within sweep, in matrix zavp.realStim) or number of trough in brst matrix

segms = zeros(numel(vertcat(zavp.realStim(:).r)), 4);%preallocation of memory for stimuli moments
z = 1;%through counter of stimuli
for sw = 1:hd.lActualEpisodes
    for t = 1:numel(zavp.realStim(sw).r) %run over synchro-events
        segms(z, 1) = zavp.realStim(sw).r(t) / zavp.rarStep;%position of synchro-event (ms from record begin (from 0))
        %segms(z, 1) = segms(z, 1) + sepOnsetPeak(15, sw).r(z, 1); disp(' shifted start point ') %shift to SEP onset
        segms(z, 2) = zavp.stimCh;%number of channel where syncro-event was detected
        segms(z, 3) = sw;%number of sweep where syncro-event was detected
        segms(z, 4) = t;%number of stimulus (in matrix zavp.realStim) or number of trough in brst matrix
        z = z + 1;%through counter of stimuli
    end
end
segms(z:end, :) = [];%delete excess

tW = round(size(lfp, 1) ./ size(segms, 1));%length of view window
if isempty(segms) %isempty(zavp.realStim(1).r) %no stimulation
    tW = 10e3;%time window
    segms = zeros(floor(size(lfp, 1) / tW), 4);%10s segments
    z = 1;
    segms(z, 1) = 0;%position of synchro-event (ms)
    segms(z, 2) = NaN;%number of channel where syncro-event was detected
    segms(z, 3) = 1;%number of sweep where syncro-event was detected
    segms(z, 4) = -1;%number of stimulus (in matrix zavp.realStim) or number of trough in brst matrix
    z = 2;
    while ((segms(z - 1, 1) + tW) < size(lfp, 1))
        segms(z, 1) = segms(z - 1, 1) + tW;%position of synchro-event (ms)
        segms(z, 2) = NaN;%number of channel where syncro-event was detected
        segms(z, 3) = 1;%number of sweep where syncro-event was detected
        segms(z, 4) = -1;%number of stimulus (in matrix zavp.realStim) or number of trough in brst matrix
        z = z + 1;
    end
    segms(z:end, :) = [];%delete excess
end

timeStr = repmat(' ', size(segms, 1), 8);%"Moscow" time of the segment begin (part of day)
for t = 1:size(segms, 1) %run over segments
    timeStr(t, :) = MoscowTime(hd.recTime(1) + (segms(t, 1) * 1e-3));%"Moscow" time of the segment begin (part of day)
end

if isappdata(evi, 'brst') %remove previous data
    rmappdata(evi, 'brst')%remove previous data
end
set(handles.ShowBrst, 'Enable', 'off')
if exist([pth, flNm, '_brst.mat'], 'file')
    load([pth, flNm, '_brst.mat'])
    if ~isstruct(brst)
        jj = brst;
        clear brst
        brst(1:size(jj, 1)) = struct('t', [], 'ch', []);
        for t = 1:size(jj, 1)
            brst(t).t = jj(t, 1:2);%begin-end time moments of bursts
            brst(t).ch = 1:hd.nADCNumChannels;%channels of bursts
        end
    end
    setappdata(evi, 'brst', brst);
    setappdata(evi, 'selectedBrst', []);%selected burst
    set(handles.ShowBrst, 'Enable', 'on')
end


setappdata(evi, 'lfp', lfp)
setappdata(evi, 'spks', spks)
setappdata(evi, 'spksOrig', spksOrig)
setappdata(evi, 'lfpVar', lfpVar)
setappdata(evi, 'hd', hd)
setappdata(evi, 'zavp', zavp)
setappdata(evi, 'segms', segms)
setappdata(evi, 'timeStr', timeStr)
set(handles.CurrSegm, 'String', '1')
set(handles.SegmTxt, 'String', num2str(size(segms, 1)))

% z = get(handles.FrstCh, 'Value');%first channel
% t = get(handles.LastCh, 'Value');%last channel
% set(handles.FrstCh, 'Value', 1);
% set(handles.LastCh, 'Value', 1);
% set(handles.FrstCh, 'String', (1:hd.nADCNumChannels)');
% set(handles.LastCh, 'String', (1:hd.nADCNumChannels)');
% if ismember(z, (1:hd.nADCNumChannels))
%     set(handles.FrstCh, 'Value', z);
% else
%     z = 1;
% end
% if ismember(t, (1:hd.nADCNumChannels))
%     set(handles.LastCh, 'Value', t);
% else
% 	set(handles.LastCh, 'Value', max(z, hd.nADCNumChannels - 5));
% end
tabData = get(handles.ShowChannels, 'Data');
if isempty(horzcat(tabData{:, 1})) %first initiation
    tabData = cell(hd.nADCNumChannels, 2);
    for t = 1:hd.nADCNumChannels %run over channels
        tabData{t, 2} = false;%channel deselected
    end
else %get previous settings
    if (size(tabData, 1) >= hd.nADCNumChannels) 
        tabData = tabData(1:hd.nADCNumChannels, :);
    else
        z = size(tabData, 1);
        for t = (z + 1):hd.nADCNumChannels %run over channels
            tabData{t, 2} = false;%channel deselected
        end
    end
end
if ~isfield(hd, 'recChNames')
    hd.recChNames = cell(hd.nADCNumChannels, 1);
    for t = 1:hd.nADCNumChannels %run over channels
        tabData{t, 1} = num2str(t);%name of channel
    end
else
    for t = 1:hd.nADCNumChannels %run over channels
        tabData{t, 1} = hd.recChNames{t};%name of channel
    end
end
if ~any(vertcat(tabData{:, 2}))
    tabData{1, 2} = true;%select only first channel
end
set(handles.ShowChannels, 'Data', tabData)%set table data
if (isequal(get(handles.BefStim, 'String'), '0') && isequal(get(handles.AftStim, 'String'), '0'))
    set(handles.BefStim, 'String', '0')
    set(handles.AftStim, 'String', num2str(tW))%length of view window
end
PlotAll(handles)


function BefStim_Callback(hObject, eventdata, handles)
% hObject    handle to BefStim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

PlotAll(handles);


function AftStim_Callback(hObject, eventdata, handles)
% hObject    handle to AftStim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

PlotAll(handles);

% --- Executes on selection change in FrstCh.
function FrstCh_Callback(hObject, eventdata, handles)
% hObject    handle to FrstCh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

PlotAll(handles);

% --- Executes on selection change in LastCh.
function LastCh_Callback(hObject, eventdata, handles)
% hObject    handle to LastCh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

PlotAll(handles);


% --- Executes on button press in NSegm.
function PNSegm_Callback(hObject, eventdata, handles)
% hObject    handle to NSegm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

evi = handles.Eview;
crSeg = handles.CurrSegm;
segms = getappdata(evi, 'segms');%number of signal segments
cs = str2double(get(crSeg, 'String'));%current segment
if isequal(hObject, handles.NSegm)
    cs = cs + str2num(get(handles.Step, 'String'));
elseif isequal(hObject, handles.PSegm)
    cs = cs - str2num(get(handles.Step, 'String'));
end
if (cs < 1)
    cs = 1;
elseif (cs > size(segms, 1))
    cs = size(segms, 1);
end
set(crSeg, 'String', num2str(cs))
PlotAll(handles);



function CurrSegm_Callback(hObject, eventdata, handles)
% hObject    handle to CurrSegm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

PlotAll(handles);



function PlotAll(handles)
evi = handles.Eview;%gui handle
lfp_ax = handles.LFP_ax;%LFP axis handle
spk_ax = handles.Spk_ax;%spikes axis handle

dataCurs = getappdata(handles.Eview, 'dataCurs');%data cursors
for t = 1:length(dataCurs)
    delete(dataCurs(t).t);%delete data cursor
    delete(dataCurs(t).l);%delete line
end
setappdata(handles.Eview, 'dataCurs', []);%data cursors

lfp = getappdata(evi, 'lfp');
spks = getappdata(evi, 'spks');
lfpVar = getappdata(evi, 'lfpVar');
hd = getappdata(evi, 'hd');
zavp = getappdata(evi, 'zavp');
segms = getappdata(evi, 'segms');
timeStr = getappdata(evi, 'timeStr');
if isappdata(evi, 'brst')
    brst = getappdata(evi, 'brst');
end

if isfield(zavp, 'chOrd')
    lfp = lfp(:, zavp.chOrd);
    spks = spks(zavp.chOrd, :);
    lfpVar = lfpVar(zavp.chOrd, 1);%for daq-files
    %lfpVar = lfpVar(zavp.chOrd);%otherwise
end
rawData = false;%if raw data needed then rawData = 1(true)

segmEdge = [-str2num(get(handles.BefStim, 'String')), ... left shifts from synchro-point (ms)
             str2num(get(handles.AftStim, 'String'))];%right shifts from synchro-point (ms)
%rCh = get(handles.FrstCh, 'Value'):get(handles.LastCh, 'Value');%channels to be read and draw
schTab = get(handles.ShowChannels, 'Data');
rCh = 1:hd.nADCNumChannels;
rCh = rCh(horzcat(schTab{:, 2}));

if (get(handles.Average, 'Value') || get(handles.WinAvr, 'Value')) %average segments
    z = str2double(get(handles.CurrSegm, 'String'));%first segment to be averaged
    t = str2double(get(handles.AvrLen, 'String'));%averaging length (number of segments to be averaged)
    sn = max(1, z):min((z + t - 1), size(segms, 1));%selected segments
    %sn = 1:size(segms, 1);%all segments
else %single segment
    sn = str2double(get(handles.CurrSegm, 'String'));%wanted sweeps or segments
    if (sn < 1)
        sn = 1;
    elseif (sn > size(segms, 1))
        sn = size(segms, 1);
    end
end
set(handles.CurrSegm, 'String', num2str(sn(1)));
set(handles.Average, 'Value', 0);%next plot single

whtLFP = ZavSynchLFP(zavp, hd, segms(sn, :), segmEdge, lfp, rCh, rawData, 0);%lfp phased with respect to stimulus moments
% whtLFP = ZavSynchLFP(zavp, hd, segms(sn, :), segmEdge, ZavFilter(lfp, prm.Fs, 'high', 5, 2), rCh, rawData, 0);%lfp phased with respect to stimulus moments
whtLFP = mean(whtLFP, 3);%average by sweeps (segments)

if get(handles.DCshift, 'Value') %DC filtration
    ajj = mean(whtLFP(1:min(5, size(whtLFP, 1)), :), 1); whtLFP = whtLFP - repmat(ajj, size(whtLFP, 1), 1);
end
% ajj = whtLFP(:, 2);
% for t = 1:2
%     whtLFP(:, t) = whtLFP(:, t) - min(whtLFP(:, t));
%     if (abs(max(whtLFP(:, t))) > 0.1)
%         whtLFP(:, t) = 9e2 * whtLFP(:, t) / max(whtLFP(:, t));
%     end
% end
% whtLFP(:, 1) = whtLFP(:, 2) - whtLFP(:, 1);
% whtLFP(:, 2) = ajj;

if rawData, prm.Fs = 1 / zavp.siS; ...%sampling frequency for raw data
else prm.Fs = 1 / (zavp.siS * zavp.rarStep); end%sampling frequency for rare data
k = 1 + ((zavp.rarStep - 1) * double(rawData == true));%multiplier for time matrix
tm = ((segmEdge(1) * k):(segmEdge(2) * k)) / k;%time matrix for raw or downsampled data (ms)


% cla(fig2)
% plot(fig2, lfp(:, 23), lfp(:, 25), '.'), hold(fig2, 'on')
% [~, z] = min(abs(tm));
% plot(fig2, whtLFP(:, 1), whtLFP(:, 3), 'r.')
% whtLFP2 = whtLFP(:, [1, 3]);

%special channels
n = 0;%number of spec.chan.
for t = 1:n
    whtLFP(:, t) = (whtLFP(:, t) - mean(whtLFP(:, t)));%amplification of signal (piezo)
    whtLFP(:, t) = whtLFP(:, t) * max(max(abs(lfp(:, rCh((n + 1):end))))) / (2 * max(abs(lfp(:, t))));
%     whtLFP(:, t) = ZavFilter(whtLFP(:, t), prm.Fs, 'high', 2, 2);%
end

%filtration
for ch = (n + 1):size(whtLFP, 2)
    %1
%     whtLFP(:, ch) = locdetrend(whtLFP(:, ch), 1e3, [0.2, 0.05]);
    %2
%     whtLFP(:, ch) = ZavFilter(whtLFP(:, ch), prm.Fs, 'high', 1, 2);%
    %3
%     whtLFP(:, ch) = ZavSpctSubtr(whtLFP(:, ch), []);%denoise by spectral subtraction
end

if (diff(segmEdge) > 300e3) %minutes
    dCoeff(4) = 60e3;%time coefficient
    timeUnit = 'min';%units
elseif ((diff(segmEdge) > 5e3) && (diff(segmEdge) <= 300e3)) %seconds
    dCoeff(4) = 1e3;%time coefficient
    timeUnit = 's';%units
else %milliseconds
    dCoeff(4) = 1;%time coefficient
    timeUnit = 'ms';%units
end
rarStep = 1 + floor(size(whtLFP, 1) / 10e3);%plot with step
if (rarStep > 1)
    ajj = whtLFP;%copy of LFP segment
    whtLFP = zeros(length(resample(ajj(:, 1), 1, rarStep)), size(ajj, 2));
    for t = 1:size(whtLFP, 2) %run over channels
        whtLFP(:, t) = resample(ajj(:, t), 1, rarStep);%downsampled LFP
    end
    tm = interp1(1:length(tm), tm, linspace(1, length(tm), size(whtLFP, 1)));%downsampled time
end

cla(lfp_ax)
setappdata(evi, 'dFrmH', []);%drag frame handle
setappdata(evi, 'bgnXY', []);%begin point of rectangle draw
setappdata(evi, 'endXY', []);%end point of rectangle draw

pltCSD = (get(handles.PlotCSD, 'Value') > 0) && ... user choosed plot CSD
         (size(whtLFP, 1) < 1e5);%short segment
if pltCSD
    if (size(whtLFP, 1) > 1e3)
%         ajj = ZavFilter(whtLFP, 1e3, 'high', 5, 2);%remove DC
        ajj = ZavRCfilt(whtLFP, 1, 1e3, 'high');%remove DC
    else
        ajj = whtLFP - repmat(mean(whtLFP(1:min(10, size(whtLFP, 1)), :), 1), numel(tm), 1);%remove DC
    end
%     ajj = ajj ./ repmat(lfpVar(rCh, 1)' ./ min(lfpVar(rCh, 1)), numel(tm), 1);%amplitude correction

    lfpCSD = -diff(ajj, 2, 2);%CSD
    k = 3;%round(diff(segmEdge) / 100) + 1;%interval of smoothing
    for t = 1:size(lfpCSD, 2) %run over channels
        lfpCSD(:, t) = smooth(lfpCSD(:, t), k);%smoothing along time
    end
    %interpolation (along channels)
    ajj = lfpCSD;
    lfpCSD = zeros(size(ajj, 1), numel(1:0.2:(numel(rCh) - 2)));
    for t = 1:size(lfpCSD, 1) %run over time
        lfpCSD(t, :) = interp1(1:(numel(rCh) - 2), ajj(t, :), 1:0.2:(numel(rCh) - 2));%interpolation along channels
    end
    
    subplot(lfp_ax);
    imagesc(tm / dCoeff(4), 2:(numel(rCh) - 1), lfpCSD'); set(gca, 'YDir', 'reverse')%normal
    %colormap('default'); hclrb = colorbar;%colobar handle
    %set(get(hclrb, 'YLabel'), 'String', 'sink - sourse')

    caxis(str2double(get(handles.CSDmapLim, 'String')) * [-1, 1])%colorbar range
    %caxis(max(max(abs(lfpCSD))) * [-1, 1] / 2);%colorbar range
end

%plot spikes and lfp
if (numel(sn) == 1) %single sweep (segment)
    jj = segms(sn, 1) + segmEdge;%absolute time of the segment boundary ( ms, e.i. from begin of file)
    sw = segms(sn, 3);%sweep
    for t = 1:numel(rCh) %run over wanted channels
        ch = rCh(t);%scanning channel
        spkDens = spks(ch, sw).tStamp((spks(ch, sw).tStamp > jj(1)) & (spks(ch, sw).tStamp < jj(end))) - segms(sn, 1);%unites in the segment
        plot(lfp_ax, [spkDens'; spkDens'] / dCoeff(4), t - repmat([0.4; 0.1], 1, numel(spkDens)), 'r', 'LineWidth', 1.5);
    end
end
p1 = str2num(get(handles.LFPScale, 'String'));
if ~isempty(p1)
    if ~isfinite(p1)
        p1 = max(max(abs(whtLFP)));
        set(handles.LFPScale, 'String', num2str(p1))
    end
end
dCoeff(1) = -1.0; dCoeff(2) = p1;%characteristic sizes (amplitude coefficient)
dCoeff(3) = dCoeff(1) * dCoeff(2);%scaling coefficient
linHandl = plot(lfp_ax, tm / dCoeff(4), (whtLFP / dCoeff(3)) + repmat(1:numel(rCh), numel(tm), 1), 'k', 'LineWidth', 1);
set(linHandl, 'ButtonDownFcn', {@Line_ButtonDownFcn, handles})%reaction on mouse-click
for t = 1:length(linHandl)
    set(linHandl(t), 'UserData', {dCoeff, t});%save amplitude coefficient and offset
end

%bursts
stt = get(handles.ShowBrst, 'State');%plot bursts?
if ((numel(sn) == 1) && isappdata(evi, 'brst') && isequal(stt, 'on')) %single segment
    jj = vertcat(brst(:).t);
    p1 = find((jj(:, 1) < (segms(sn, 1) + segmEdge(2))) & (jj(:, 2) > (segms(sn, 1) + segmEdge(1))));
    for t = p1'
        [~, ii] = intersect(rCh, brst(t).ch);%common channels
        if ~isempty(ii) %some channels contain the burst
            jj = ((tm >= (brst(t).t(1) - segms(sn, 1))) & (tm <= (brst(t).t(2) - segms(sn, 1))));
            plot(lfp_ax, tm(jj) / dCoeff(4), (whtLFP(jj, ii) / dCoeff(3)) + repmat((rCh(ii)) - (rCh(1) - 1), sum(jj), 1), 'r', 'LineWidth', 2)
        end
    end
end

set(lfp_ax, 'YDir', 'reverse', 'XLim', segmEdge / dCoeff(4), 'YLim', [-0.05, numel(rCh) + 1.05])
set(get(lfp_ax, 'XLabel'), 'String', ['time, ', timeUnit])

chnlStr = num2str(rCh');%1) common label for channels
% chnlStr = char;%2) label for man EEG
% for t = 1:numel(rCh)
% %     %fN = regexprep(hd.recChNames{rCh(t)}, 'EEG', '');
% %     %chnlStr(t, 1:length(fN)) = fN;
%     fN = regexprep(hd.label{rCh(t)}, 'EEG', '');
%     chnlStr(t, 1:length(fN)) = fN;
% end
% z = segms(sn(end), 1) / (1e3 * 60);%absolute time (minutes)
% disp([num2str(floor(z)), 'min ', num2str(round((z - floor(z)) * 60)), 'sec'])%absolute time (minutes)

set(lfp_ax, 'YTick', 1:numel(rCh), 'YTickLabel', chnlStr)
set(get(lfp_ax, 'YLabel'), 'String', 'channels', 'Units', 'normalized', 'Position', [-0.06, 0.6, 1])
plot(lfp_ax, [0, 0], get(lfp_ax, 'YLim'), 'r', 'LineWidth', 2)%nul-time line

[~, p1] = ds2nfu([0, 0], [0, abs(1 / dCoeff(1))]);%diff(p1) = length of scale bar in normalized units
set(handles.LFPscb, 'Position', [0.025, 0.15 0, diff(p1)], 'LineWidth', 5, 'Units', 'normalized', 'TextRotation', 90, 'HeadStyle', 'none',...
    'String', [num2str(round(dCoeff(2))), ' \muV'], 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');

if (numel(sn) == 1)
    h = text(0, 0, timeStr(sn, :));%time of the current segment
    set(h, 'Parent', lfp_ax);
end

% ================================ %
% ===== spikes density + lfp ===== %
cla(spk_ax)
if get(handles.SpkShow, 'Value') %show spikes
    binStep = [str2num(get(handles.SpkBin, 'String')), str2num(get(handles.HistPrecision, 'String'))];%[bin size (ms), precision (must be integer)]
    binStep(2) = round(binStep(2));%
    bins = segmEdge(1):binStep(1):segmEdge(2);%base bins
    spkDens = zeros(numel(bins) * binStep(2), numel(rCh), numel(sn));%density of spikes in time (histogram)

    for z = 1:numel(sn) %run over requested segments
        jj = segms(sn(z), 1) + segmEdge;%time segment boundary (ms from begin of file)
        sw = segms(sn(z), 3);%sweep
        bins = (jj(1):binStep(1):jj(2));%binning of time segment
        for ch = 1:numel(rCh) %run over channels
            %fast overlaped histogram
            if ~isempty(spks(rCh(ch), sw).tStamp) %spikes detected
                for k = 1:binStep(2) %run over precision steps
                    spkDens(k:binStep(2):end, ch, z) = histc(spks(rCh(ch), sw).tStamp, bins + (k - 1) * (binStep(1) / binStep(2)));
                end
            end
        end
    end
    spkDens = mean(spkDens, 3);%mean spikes density
    bins = (segmEdge(1) + (binStep(1) / 2)) + (0:(size(spkDens, 1) - 1)) * (binStep(1) / binStep(2));%
    jj = (bins >= segmEdge(1)) & (bins <= segmEdge(2));%values within requested range only
    spkDens = spkDens(jj, :);%values within requested range only
    bins = bins(jj);%values within requested range only

    p1 = str2num(get(handles.MUAScale, 'String'));
    if ~isempty(p1)
        if ~isfinite(p1)
            p1 = max(max(spkDens));
            set(handles.LFPScale, 'String', num2str(p1))
        end
    end

    dCoeff(1) = -1; dCoeff(2) = p1;%characteristic sizes
    dCoeff(3) = dCoeff(1) * dCoeff(2);%scaling coefficient

    rawData = get(handles.SpkMap, 'Value');%plot spike density curves
    clr = 'b';
    if ((numel(rCh) > 1) && rawData) %show many channels
    %     %interpolation (smooth along channels)
    %     spkDensSC = zeros(size(spkDens, 1), numel(1:0.25:numel(rCh)));
    %     for t = 1:size(spkDens, 1) %run over time
    %         spkDensSC(t, :) = interp1(1:numel(rCh), spkDens(t, :), 1:0.25:numel(rCh));%interpolation along channels
    %     end
        spkDensSC = spkDens;
        %smooth along time
        for t = 1:size(spkDens, 2) %rum over channels
            spkDensSC(:, t) = smooth(spkDens(:, t), 3);%smooth along time
        end

        subplot(spk_ax)
        imagesc(bins / dCoeff(4), 1:numel(rCh), spkDensSC'); set(gca, 'YDir', 'reverse')%normal
        %clrMap = colormap('gray'); colormap(clrMap(end:-1:1, :))
        shading flat;
        caxis([str2double(get(handles.MUAmapMin, 'String')), str2double(get(handles.MUAmapMax, 'String'))])%colorbar range

        clr = 'w';
    end
    linHandl = plot(spk_ax, bins / dCoeff(4), (spkDens / dCoeff(3)) + repmat(1:numel(rCh), numel(bins), 1), clr, 'LineWidth', 2);
    set(spk_ax, 'YDir', 'reverse')
    set(spk_ax, 'YTick', 1:numel(rCh), 'YTickLabel', num2str(rCh'), 'YLim', [0, numel(rCh) + 1])

    set(linHandl, 'ButtonDownFcn', {@Line_ButtonDownFcn, handles})%reaction on mouse-click
    for t = 1:length(linHandl)
        set(linHandl(t), 'UserData', {dCoeff, t});%save amplitude coefficient and offset
    end

    set(get(spk_ax, 'YLabel'), 'String', 'channels')
    [~, p1] = ds2nfu([0, 0], [0, abs(1 / dCoeff(1))]);%diff(p1) = length of scale bar in normalized units
    set(handles.MUAscb, 'Position', [0.525, 0.15 0, diff(p1)], 'LineWidth', 5, 'Units', 'normalized', 'TextRotation', 90, 'HeadStyle', 'none',...
        'String', [num2str(round(dCoeff(2) * 10) / 10), 'units/', num2str(binStep(1)), 'ms'], 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');
    plot(spk_ax, [0, 0], get(gca, 'YLim'), 'r', 'LineWidth', 2)%nul-time line

    set(get(spk_ax, 'XLabel'), 'String', ['time, ', timeUnit])
    set(spk_ax, 'XLim', segmEdge / dCoeff(4))

    if (numel(sn) == 1)
        h = text(0, 0, timeStr(sn, :));%time of the current segment
        set(h, 'Parent', spk_ax);
    end
end

% bufStr = '_segment';%litter for define number of sweep or segment
% if (numel(sn) > 1)
%     bufStr = [bufStr, num2str(sn(1)), '-', num2str(sn(end))];
% else
%     bufStr = [bufStr, num2str(sn)];
% end
% p1 = find(zavp.file(1:end - 1) == '\', 1, 'last') + 1;%
% p2 = find(zavp.file == '.', 1, 'last') - 1; if (isempty(p2) || (p1 >= p2)), p2 = length(zavp.file) - 1; end
% fN = [zavp.file(p1:p2), bufStr];
% set(evi, 'Name', ['Eview - ', fN])
% set(gcf, 'FileName', ['c:\Users\Zav\Pictures\', fN, '.emf'], 'Name', fN, 'NumberTitle', 'off')
% fN = regexprep(fN, '_', '\\_');%replace character '_' by '\_'
% title(fN)



% --- Executes on button press in PlotCSD.
function PlotCSD_Callback(hObject, eventdata, handles)
% hObject    handle to PlotCSD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(hObject, 'Value')
    set(handles.CSDmapLim, 'Enable', 'on')
else
    set(handles.CSDmapLim, 'Enable', 'off')
end
PlotAll(handles)


% --- Executes on button press in Average.
function Average_Callback(hObject, eventdata, handles)
% hObject    handle to Average (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(hObject, 'Value', 1);
PlotAll(handles)


% --- Executes on button press in SpkMap.
function SpkMap_Callback(hObject, eventdata, handles)
% hObject    handle to SpkMap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(hObject, 'Value')
    set(handles.MUAmapMax, 'Enable', 'on')
    set(handles.MUAmapMin, 'Enable', 'on')
else
    set(handles.MUAmapMax, 'Enable', 'off')
    set(handles.MUAmapMin, 'Enable', 'off')
end
PlotAll(handles)


function LFPScale_Callback(hObject, eventdata, handles)
% hObject    handle to LFPScale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

PlotAll(handles)


function SpkPrg_Callback(hObject, eventdata, handles)
% hObject    handle to SpkPrg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

evi = handles.Eview;
prg = str2num(get(handles.SpkPrg, 'String'));%threshold multiplier
prg = prg(1);%threshold multiplier
set(handles.SpkPrg, 'String', num2str(prg))
spksOrig = getappdata(evi, 'spksOrig');
lfpVar = getappdata(evi, 'lfpVar');
hd = getappdata(evi, 'hd');
spks = spksOrig;
if isfield(spksOrig, 'ampl')
    for sw = 1:hd.lActualEpisodes
        for ch = 1:hd.nADCNumChannels
            ii = double(spksOrig(ch, sw).ampl) <= (-lfpVar(ch) * prg);
            spks(ch, sw).tStamp = spksOrig(ch, sw).tStamp(ii);
            spks(ch, sw).ampl = spksOrig(ch, sw).ampl(ii);
        end
    end
end
setappdata(evi, 'spks', spks)
PlotAll(handles)

function timeStr = MoscowTime(timeNum)
%convert time as number to string in format hh:mm
%INPUTS
%timeNum - time as number (seconds from day beginning)
%OUTPUTS
%timeStr - time as string in format hh:mm

tm = zeros(1, 2);
tm(1) = floor(timeNum / 3600);%hours
tm(2) = floor((timeNum - (tm(1) * 3600)) / 60);%minutes
tm(3) = floor(timeNum - (tm(1) * 3600) - tm(2) * 60);%seconds

timeStr = num2str(tm(1));
if (tm(1) < 10)
    timeStr = ['0', timeStr];
end
timeStr = [timeStr, ':'];
if (tm(2) < 10)
    timeStr = [timeStr, '0'];
end
timeStr = [timeStr, num2str(tm(2))];
timeStr = [timeStr, ':'];
if (tm(3) < 10)
    timeStr = [timeStr, '0'];
end
timeStr = [timeStr, num2str(tm(3))];


function SpkBin_Callback(hObject, eventdata, handles)
% hObject    handle to SpkBin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

PlotAll(handles)


% --- Executes on mouse press over axes background.
function Line_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to LFP_ax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

evi = handles.Eview;
dataCurs = getappdata(evi, 'dataCurs');%data cursors
segms = getappdata(evi, 'segms');%segments of signal
timeStr = getappdata(evi, 'timeStr');%synchromoments
if isempty(dataCurs)
    dataCurs = struct('t', [], 'l', []);
    nCrs = 1;%number of data cursors
else
    nCrs = length(dataCurs) + 1;%number of data cursors
end

%creat new data cursor
xy = get(get(hObject, 'Parent'), 'CurrentPoint');
xy = xy(1, 1:2);
[p1, p2] = ds2nfu(xy(1), xy(2));
dCoeff = get(hObject, 'UserData');%amplitude coefficient
xy(2) = dCoeff{1}(3) * (xy(2) - dCoeff{2});

if isequal(get(get(handles.LFP_ax, 'XLabel'), 'String'), 'time, s')
    dCoeff = 1e3;
else
    dCoeff = 1;
end
cs = str2double(get(handles.CurrSegm, 'String'));%current segment
[~, m] = min(abs(segms(:, 1) - (segms(cs, 1) + xy(1) * dCoeff)));%number of the nearest timestamp

dataCurs(nCrs).t = annotation('textbox', 'Units', 'normalized', 'Position', [p1, 1 - p2 - 0.01, 0.055, 0.045], ...
    'String', {['X:', num2str(round(xy(1) * 10) / 10)], ['Y:', num2str(round(xy(2) * 10) / 10)], timeStr(m, :)}, ...
    'BackgroundColor', [1, 1, 1], 'EdgeColor', [0, 1, 0], 'UserData', {xy, hObject}, 'FontSize', 8);
dataCurs(nCrs).l = annotation('line', 'Units', 'normalized', 'Position', [p1, 1 - p2 - 0.01, 0.0, -0.025], ...
    'Color', [0, 1, 0]);

if (nCrs > 1)
    psn = zeros(nCrs, 3);
    for t = 1:nCrs
        userdata = get(dataCurs(t).t, 'UserData');
        psn(t, 1:2) = userdata{1};%clicked point (data space)
        psn(t, 3) = userdata{2};%clicked line (handle)
    end
    hstH = unique(psn(:, 3));
    for k = 1:numel(hstH)
        jj = (psn(:, 3) == hstH(k));
        psnL = psn(jj, 1:2);
        hL = vertcat(dataCurs(jj).t);
        [~, ii] = sort(psnL(:, 1));
        for z = numel(ii):-1:2
            dx = round((psnL(ii(z), 1) - psnL(ii(z - 1), 1)) * 10) / 10;
            dy = round((psnL(ii(z), 2) - psnL(ii(z - 1), 2)) * 10) / 10;
            [~, m] = min(abs(segms(:, 1) - (segms(cs, 1) + xy(1) * dCoeff)));%number of the nearest timestamp
            set(hL(ii(z)), 'String', {['dX:', num2str(dx)], ...
                                      ['dY:', num2str(dy)], ...
                                      timeStr(m, :)});
        end
    end
end
setappdata(evi, 'dataCurs', dataCurs);

% --------------------------------------------------------------------
function ShowBrst_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to ShowBrst (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

PlotAll(handles);


% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function Eview_WindowButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to Eview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

evi = handles.Eview;%main form handle
lfpAx = handles.LFP_ax;%LFP axes
% spkAx = handles.Spk_ax;%spikes axes

x_lim = get(lfpAx, 'XLim');
y_lim = get(lfpAx, 'YLim');
bgnXY = get(lfpAx, 'CurrentPoint');
clkIn = zeros(1, 3);
clkIn(1) = ((bgnXY(1, 1) >= x_lim(1)) && (bgnXY(1, 1) <= x_lim(2)) && ...
            (bgnXY(1, 2) >= y_lim(1)) && (bgnXY(1, 2) <= y_lim(2))); %click point is in axes limits

% x_lim = get(spkAx, 'XLim');
% y_lim = get(spkAx, 'YLim');
% bgnXY = get(spkAx, 'CurrentPoint');
% clkIn(2) = ((bgnXY(1, 1) >= x_lim(1)) && (bgnXY(1, 1) <= x_lim(2)) && ...
%             (bgnXY(1, 2) >= y_lim(1)) && (bgnXY(1, 2) <= y_lim(2))); %click point is in axes limits
        
if ~any(clkIn) %no one axes clicked
    setappdata(evi, 'bgnXY', []);%start point of frame
    return
end
if clkIn(1) %click on signalAx
    axN = lfpAx;%clicked axes
% elseif clkIn(2) %click on parametrAx
%     axN = spkAx;%clicked axes
else
    axN = [];%error
end
dFrmH = getappdata(evi, 'dFrmH');%handle of frame
axFrm = getappdata(evi, 'axFrm');%handle of clicked axes
bgnXY = get(axN, 'CurrentPoint');
setappdata(evi, 'bgnXY', bgnXY);%start point of frame
setappdata(evi, 'endXY', []);%end point of frame
if (isempty(dFrmH) || ~isequal(axFrm, axN)) %new frame (may on other axis)
    if (~isempty(dFrmH) && ishandle(dFrmH)) %frame exist on other axes
        delete(dFrmH)%destroy object
    end
    x_lim = get(axN, 'XLim');
    y_lim = get(axN, 'YLim');
    dFrmH = plot(axN, mean(x_lim), mean(y_lim), 'r-');
    set(axN, 'XLim', x_lim, 'YLim', y_lim);
    setappdata(evi, 'dFrmH', dFrmH);%handle of frame
    setappdata(evi, 'axFrm', axN);%handle of clicked axes
end


% --- Executes on mouse motion over figure - except title and menu.
function Eview_WindowButtonMotionFcn(hObject, eventdata, handles)
% hObject    handle to Eview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

evi = handles.Eview;
bgnXY = getappdata(evi, 'bgnXY');%start point of frame
endXY = getappdata(evi, 'endXY');%end point of frame
if (~isempty(bgnXY) && isempty(endXY))
    axFrm = getappdata(evi, 'axFrm');%handle of axes with frame
    endXY = get(axFrm, 'CurrentPoint');
    x_lim = get(axFrm, 'XLim');
    y_lim = get(axFrm, 'YLim');
    if ((endXY(1, 1) >= x_lim(1)) && (endXY(1, 1) <= x_lim(2)) && ...
        (endXY(1, 2) >= y_lim(1)) && (endXY(1, 2) <= y_lim(2))) %click point is in axes limits
        dFrmH = getappdata(evi, 'dFrmH');%drag frame handle
        set(dFrmH, 'Xdata', [bgnXY(1, 1), endXY(1, 1), endXY(1, 1), bgnXY(1, 1), bgnXY(1, 1)], ...
            'Ydata', [bgnXY(1, 2), bgnXY(1, 2), endXY(1, 2), endXY(1, 2), bgnXY(1, 2)])
    end
end


% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function Eview_WindowButtonUpFcn(hObject, eventdata, handles)
% hObject    handle to Eview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

evi = handles.Eview;
spks = getappdata(evi, 'spks');
bgnXY = getappdata(evi, 'bgnXY');%left upper point
if (~isempty(spks) && ~isempty(bgnXY))
    axFrm = getappdata(evi, 'axFrm');%handle of axes with frame
    endXY = get(axFrm, 'CurrentPoint');%right lower point
%     setappdata(evi, 'endXY', endXY);
%     dFrmH = getappdata(evi, 'dFrmH');%drag frame handle
%     set(dFrmH, 'Xdata', [bgnXY(1, 1), endXY(1, 1), endXY(1, 1), bgnXY(1, 1), bgnXY(1, 1)], ...
%         'Ydata', [bgnXY(1, 2), bgnXY(1, 2), endXY(1, 2), endXY(1, 2), bgnXY(1, 2)])
    x_lim = sort([bgnXY(1, 1), endXY(1, 1)]);%horizontal range
    y_lim = sort([bgnXY(1, 2), endXY(1, 2)]);%vertical range
    hd = getappdata(evi, 'hd');%header
    segms = getappdata(evi, 'segms');%segments
    sn = str2num(get(handles.CurrSegm, 'String'));%segments of signal (sweeps)
    
    %rCh = get(handles.FrstCh, 'Value'):get(handles.LastCh, 'Value');%channels to be read and draw
    schTab = get(handles.ShowChannels, 'Data');
    rCh = 1:hd.nADCNumChannels;
    rCh = rCh(horzcat(schTab{:, 2}));%selected channels
    
    lfpChlds = get(handles.LFP_ax, 'Children');
    z = 3;
    dCoeff = [];
    while isempty(dCoeff)
        dCoeff = get(lfpChlds(z), 'UserData');
        z = z + 1;
    end
    dCoeff = dCoeff{1};%scale factors
    
    bT.t = segms(sn, 1) + x_lim * dCoeff(4);%begin-end of currently defined burst (or highlighted area)
    jj = ceil(y_lim(1)):floor(y_lim(2));%highlighted vertical range
    jj = jj((jj > 0) & (jj <= numel(rCh)));%selected channels
    bT.ch = rCh(jj);%selected channels
    bT.t(1) = max([1, bT.t(1)]);
    bT.t(2) = min([bT.t(2), (diff(hd.recTime) * 1e3) - 1]);
    if isequal(get(handles.PwrSpectrWin, 'State'), 'on') %spectrum in window selected
        %immediately remove the frame
        delete(getappdata(evi, 'dFrmH'))
        setappdata(evi, 'dFrmH', []);%drag frame handle
        setappdata(evi, 'bgnXY', []);%begin point of rectangle draw
        setappdata(evi, 'endXY', []);%end point of rectangle draw
        
        %plot spectra
        if (~isempty(bT.ch) && (diff(bT.t) > 5)) %traces selected
            figure
            lfp = getappdata(evi, 'lfp');
            zavp = getappdata(evi, 'zavp');
            hd = getappdata(evi, 'hd');
            whtLFP = ZavSynchLFP(zavp, hd, [bT.t(1), NaN, 1, 1], [0, round(diff(bT.t))], lfp, bT.ch, 0, 0);
            whtLFP = ZavFilter(whtLFP, 1e3, 'high', 5, 2);

            if isfield(zavp, 'dwnSmplFrq') %field exist
                prm.Fs = zavp.dwnSmplFrq;%discretization frequency for signals to be treated (Hz)
            else
                prm.Fs = 1000;%discretization frequency for signals to be treated (Hz)
            end
            prm.fpass = [5, 100];%frequency range (predominantly for sensory response)
            prm.tapers = [3, 5];%tapers
            prm.pad = 2;%padding factor
            prm.trialave = 0;%channel by channel, sweep by sweep (don't average when 0)
            
            [pwrSpct, frq] = mtspectrumc(whtLFP, prm);%mean power spectrum (mean by sweeps)
            hold on
            plot(frq, pwrSpct)
            plot(frq, mean(pwrSpct, 2), 'k', 'LineWidth', 2)
            xlim(frq([1, end]))
            lgndStr = {};
            for z = bT.ch
                lgndStr{end + 1} = ['ch', num2str(z)];
            end
            lgndStr{end + 1} = 'average';
            legend(lgndStr)
        end
    else %complete selected burst
        if isappdata(evi, 'brst')
            brst = getappdata(evi, 'brst');
            z = length(brst) + 1;
        else
            brst = struct('t', [], 'ch', []);
            z = 1;
        end
    %     if isequal(axFrm, lfpAx) %LFP axes
    %     elseif isequal(axFrm, spkAx) %frame is on parameters axes
    %     end
        %immediately remove the frame
        delete(getappdata(evi, 'dFrmH'))
        setappdata(evi, 'dFrmH', []);%drag frame handle
        setappdata(evi, 'bgnXY', []);%begin point of rectangle draw
        setappdata(evi, 'endXY', []);%end point of rectangle draw

        sb = false;%if any burst exist on selected area
        b = [];%number of bursts within selected area
        jj = vertcat(brst(:).t);%all detected bursts
        if ~isempty(jj) %bursts exist
            b = find((jj(:, 1) <= bT.t(1)) & (jj(:, 2) >= bT.t(1)));%left board on burst
            b = [b, find((jj(:, 1) <= bT.t(2)) & (jj(:, 2) >= bT.t(2)))];%right board on burst
            if ~isempty(b) %burst exist within selected area
                sb = true;%some burst exist on selected range
            else %bouth b1 and b2 are empty
                b = find((jj(:, 1) >= bT.t(1)) & (jj(:, 1) <= bT.t(2)));%number of bursts within selected area
                if ~isempty(b) %burst exist within selected area
                    sb = true;%some burst exist on selected range
                end
            end
        end
        if ~isempty(b) %burst exist within selected area
            if ~any(ismember(horzcat(brst(b).ch), bT.ch))
                sb = false;%no bursts selected
            end
        end
        if (sb) %some burst exist on selected range
            setappdata(evi, 'selectedBrst', b);%selected burst
        else %the new burst demand
            if ((diff(round(bT.t)) > 5) && (numel(bT.ch) > 0))
                brst(z).t = round(bT.t);%absolute time (ms = samples)
                brst(z).ch = bT.ch;%channels
                setappdata(evi, 'selectedBrst', z);%selected burst
                setappdata(evi, 'brst', brst);
                set(handles.ShowBrst, 'Enable', 'on')
                set(handles.ShowBrst, 'State', 'on')
                PlotAll(handles);
            end
        end
    end
end


% --------------------------------------------------------------------
function DeleteBrst_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to DeleteBrst (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

evi = handles.Eview;
if isappdata(evi, 'brst') %bursts exist
    slB = getappdata(evi, 'selectedBrst');%selected burst
    if ~isempty(slB) %any burst selected
        brst = getappdata(evi, 'brst');
        brst(slB) = [];%delete selected burst
        if isempty(brst)
            rmappdata(evi, 'brst');
            set(handles.ShowBrst, 'Enable', 'off')
        else
            setappdata(evi, 'brst', brst);
        end
        setappdata(evi, 'selectedBrst', []);
        PlotAll(handles);
    end
end


% --- Executes when user attempts to close Eview.
function Eview_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to Eview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

evi = handles.Eview;%main form handle
if isappdata(evi, 'brst')
    ynBtn = questdlg('Save resultes?', 'Save on exit', 'Yes', 'No', 'No');%ask what to do
    if isequal(ynBtn, 'Yes') %save requested
        SaveBursts_ClickedCallback([], [], handles);%save bursts
    end
end
% Hint: delete(hObject) closes the figure
delete(hObject);


% --------------------------------------------------------------------
function CSDcalcul_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to CSDcalcul (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

evi = handles.Eview;%main form handle

%calculate parameters of CSD (Cortical Spreading Depression)
if isappdata(evi, 'brst')
    pth = getappdata(evi, 'pth');%directory
    flNm = getappdata(evi, 'flNm');%filename
    
    lfp = getappdata(evi, 'lfp');
    %lfp = double(lfp);
    spks = getappdata(evi, 'spks');
    %lfpVar = getappdata(evi, 'lfpVar');
    hd = getappdata(evi, 'hd');
    zavp = getappdata(evi, 'zavp');
    brst = getappdata(evi, 'brst');
    prm(1:length(brst)) = struct('p', []);
    jj = vertcat(brst(:).t);
    [~, ii] = sort(jj(:, 1));
    brst = brst(ii);%sorted bursts
    
    filtFlag = false(hd.nADCNumChannels, 1);%true if lfp was fitered
    for nch = 1:hd.nADCNumChannels
        try
            cscHd = Nlx2MatCSC([pth, flNm, '\', hd.recChNames{nch}, '.ncs'], [0 0 0 0 0], 1, 1, []);%header (lfp)
        catch
            ynBtn = questdlg('Low-cut filter enabled?', 'low-cut filter', 'Yes', 'No', 'No');%ask what to do
            if isequal(ynBtn, 'Yes') %save requested
                filtFlag = true(hd.nADCNumChannels, 1);%all is true (filter was enabled)
            else
                filtFlag = false(hd.nADCNumChannels, 1);%all is false (filter was not enabled)
            end
            break;%out of (for nch)
        end
        if ~isempty(strfind(horzcat(cscHd{:}), '-DSPLowCutFilterEnabled True')) %lowcut filter enabeled
            filtFlag(nch) = true;%lowcut-filtered lfp
        end
    end
    for t = 1:length(brst)
        prm(t).p = Inf(length(brst(t).ch), 5);%parameters
        for nch = 1:length(brst(t).ch)
            rCh = brst(t).ch(nch);%current channel
            nt = [brst(t).t(1), brst(t).t(2)];%noses and tails
            if filtFlag(rCh) %lowcut filter enabeled
                prg = [mean(lfp(:, rCh)), std(lfp(:, rCh) / 2)];%threshold
                nose_tail = lfp(nt(1) + (0:2e3), rCh);%"nose" of CSD
                while (abs(mean(nose_tail) - prg(1)) > prg(2)) %null-line is reached
                    nt(1) = max(nt(1) - 2e3, 1);%left side
                    nose_tail = lfp(nt(1) + (0:2e3), rCh);%"nose" of CSD
                    if (nt(1) == 1) %too small value
                        break;%out of (while)
                    end
                end
                
                nose_tail = lfp(nt(2) + (-2e3:0), rCh);%"tail" of CSD
                while ((abs(mean(nose_tail) - prg(1)) > prg(2)) || (sum(lfp(nt(1):nt(2), rCh)) < prg(2))) %null-line is reached
                    nt(2) = min(nt(2) + 2e3, size(lfp, 1));%left side
                    nose_tail = lfp(nt(2) + (-2e3:0), rCh);%"tail" of CSD
                    if (nt(2) == size(lfp, 1)) %too large value
                        break;%out of (while)
                    end
                end

                ii = nt(1):nt(2);%range
                whtLFP = cumsum(double(lfp(ii, rCh)));%recovery non-filtered signal
            else
                prg = [mean(lfp(:, rCh)), std(lfp(:, rCh))];%threshold
                nose_tail = lfp(nt(1) + (0:2e3), rCh);%"nose" of CSD
                while (abs(mean(nose_tail) - prg(1)) > prg(2)) %null-line is reached
                    nt(1) = max(nt(1) - 2e3, 1);%left side
                    nose_tail = lfp(nt(1) + (0:2e3), rCh);%"nose" of CSD
                    if (nt(1) == 1) %too small value
                        break;%out of (while)
                    end
                end
                
                nose_tail = lfp(nt(2) + (-2e3:0), rCh);%"tail" of CSD
                while (abs(mean(nose_tail) - prg(1)) > prg(2)) %null-line is reached
                    nt(2) = min(nt(2) + 2e3, size(lfp, 1));%left side
                    nose_tail = lfp(nt(2) + (-2e3:0), rCh);%"tail" of CSD
                    if (nt(2) == size(lfp, 1)) %too large value
                        break;%out of (while)
                    end
                end
                
                ii = nt(1):nt(2);%range
                whtLFP = ZavFilter(double(lfp(ii, rCh)), prm.Fs, 'low', 3, 2);%smooth(lfp(jj, rCh), 31);
                whtLFP = whtLFP - mean(whtLFP(1:2e3));%adduction to zeros
            end
            [a, p] = min(whtLFP);%CSD peak (negative only)
            prm(t).p(nch, 2) = p;%time of the first (negative) peak (ms from frame begin)
            prm(t).p(nch, 3) = a;%amplitude of the first (negative) peak (uV)
            
            da = prm(t).p(nch, 3) * 0.2;%20 percent of negative amplitude
            p = find(whtLFP(1:prm(t).p(nch, 2)) > da, 1, 'last');
            if ~isempty(p)
                prm(t).p(nch, 1) = p;%onset (CSD start, ms from frame begin)
            end
            
            p = find(whtLFP(prm(t).p(nch, 2):end) < da, 1, 'last');
            if ~isempty(p)
                prm(t).p(nch, 4) = p + prm(t).p(nch, 2) - 1;%offset (CSD end, ms from frame begin)
            end
            
            da = prm(t).p(nch, 3) * 0.5;%50 percent of negative amplitude
            p = find(whtLFP(prm(t).p(nch, 2):-1:1) > da, 1, 'first');%left half high
            if ~isempty(p)
                a = find(whtLFP(prm(t).p(nch, 2):end) > da, 1, 'first');%rigth half high
                if ~isempty(a)
                    prm(t).p(nch, 5) = a + p;%width on half high
                end
            end
            %clf, hold on, jj = prm(t).p(nch, [1, 2, 4]); jj = jj(isfinite(jj)); plot(ii - ii(1) + 1, whtLFP, jj, whtLFP(jj), 'or'), jj = prm(t).p(nch, 2); if isfinite(jj), plot(jj + [-p, a], whtLFP(jj + [-p, a]), 'gs'), end
            
            prm(t).p(nch, [1, 2, 4]) = prm(t).p(nch, [1, 2, 4]) + ii(1) - 1;%time of anythink from record begin (ms)
            
            %whtLFP = cumsum(lfp(:, rCh)); clf, hold on, jj = prm(t).p(nch, [1, 2, 4]); jj = jj(isfinite(jj)); plot(1:size(lfp, 1), whtLFP, jj, whtLFP(jj), 'or'), jj = prm(t).p(nch, 2); if isfinite(jj), plot(jj + [-p, a], whtLFP(jj + [-p, a]), 'gs'), end
            %clf, hold on, jj = prm(t).p(nch, [1, 2, 4]); jj = jj(isfinite(jj)); plot(1:size(lfp, 1), lfp(:, rCh), jj, lfp(jj, rCh), 'or'), jj = prm(t).p(nch, 2); if isfinite(jj), plot(jj + [-p, a], lfp(jj + [-p, a], rCh), 'gs'), end
        end
    end

    %export parameter of CSD to file
    fid = fopen([pth, flNm, '_CSDprm.dat'], 'w');

    fprintf(fid, ['MUA threshold\t', get(handles.SpkPrg, 'String'), '\n']);%threshold multiplier
    fprintf(fid, ['Vertical scale (uV)\t', get(handles.LFPScale, 'String'), '\n']);%vertical scale
    fprintf(fid, 'Window length (ms)\t%d\n', str2num(get(handles.AftStim, 'String')) + str2num(get(handles.BefStim, 'String')));%time window length
    fprintf(fid, ['Step\t', get(handles.Step, 'String'), '\n\n']);%step (number of windows)

    spkDens = zeros(hd.nADCNumChannels, 1);%total MUA on channels
    p1 = zeros(hd.nADCNumChannels, 1);%current begin
    len = zeros(hd.nADCNumChannels, 1);%total lengthes of interburst intervals
    for t = 1:length(brst)
        for rCh = 1:hd.nADCNumChannels
            if ismember(rCh, brst(t).ch)
                p2 = brst(t).t(1);%current right edge
                spkDens(rCh) = spkDens(rCh) + sum((spks(rCh).tStamp >= p1(rCh)) & (spks(rCh).tStamp <= p2));%MUA
                len(rCh) = len(rCh) + (p2 - p1(rCh));%total length (ms)
                p1(rCh) = brst(t).t(2);%new left edge
            end
        end
    end
    p2 = size(lfp, 1);%common right edge
    for rCh = 1:hd.nADCNumChannels
        spkDens(rCh) = spkDens(rCh) + sum((spks(rCh).tStamp >= p1(rCh)) & (spks(rCh).tStamp <= p2));%MUA
        len(rCh) = len(rCh) + (p2 - p1(rCh));%total length
        spkDens(rCh) = spkDens(rCh) / (len(rCh) * 1e-3);%MUA per second
    end
    fprintf(fid, 'channels\tMUA background (1/s)\n');
    for rCh = 1:hd.nADCNumChannels
        fprintf(fid, '%d\t%3.3f\n', [rCh, spkDens(rCh)]);
    end
    fprintf(fid, '\n\n');

    fprintf(fid, 'N\tbegin\tend\tchannel\tonset\tpeak\tampl\toffset\tWidth\tMUA\n');
    fprintf(fid, '#\tms\tms\t#\tms\tms\tuV\tms\tms\tunits\n');
    for t = 1:length(brst)
        fprintf(fid, '%d\t%d\t%d\n', [t, brst(t).t]);
        for rCh = 1:hd.nADCNumChannels
            fprintf(fid, '\t\t\t%d\t', rCh);
            if ismember(rCh, brst(t).ch)
                nch = find((brst(t).ch == rCh));
                fprintf(fid, '%d\t%d\t%6.2f\t%d\t%d\t', prm(t).p(nch, :));
                spkDens = NaN;
                if all(isfinite(prm(t).p(nch, 1:2)))
                    spkDens = sum((spks(rCh).tStamp >= (prm(t).p(nch, 1) - 1e4)) & (spks(rCh).tStamp <= (prm(t).p(nch, 2) + 0)));%MUA between 10 ms befor CSD onset and CSD peak
                end
                    fprintf(fid, '%d\n', spkDens);
            else
                fprintf(fid, '\t\t\t\t\t\t\n');
            end
        end
        fprintf(fid, '\n');
    end
    fclose(fid);
end

% --------------------------------------------------------------------
function DetectBrst_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to DetectBrst (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

evi = handles.Eview;%main form handle
%choose channel


% --- Executes on button press in DCshift.
function DCshift_Callback(hObject, eventdata, handles)
% hObject    handle to DCshift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

PlotAll(handles);


% --------------------------------------------------------------------
function SaveBursts_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to SaveBursts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%save bursts

evi = handles.Eview;%main form handle
pth = getappdata(evi, 'pth');%directory
flNm = getappdata(evi, 'flNm');%filename
brst = getappdata(evi, 'brst');
psn = get(evi, 'Position');
hWt = waitbar(0.5, 'saving bursts', 'Name', 'Save', 'Units', 'pixels', 'Position', [psn(1), psn(2) + psn(4) - 75, 360, 75]);%waitbar

jj = vertcat(brst(:).t);
[~, ii] = sort(jj(:, 1));
brst = brst(ii);%sorted bursts

if exist([pth, flNm, '_brst.mat'], 'file')
    save([pth, flNm, '_brst.mat'], 'brst', '-append')
else
    save([pth, flNm, '_brst.mat'], 'brst')
end
close(hWt);%destroy waitbar



function HistPrecision_Callback(hObject, eventdata, handles)
% hObject    handle to HistPrecision (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

PlotAll(handles);


% --------------------------------------------------------------------
function ExportImages_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to ExportImages (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%export current images
evi = handles.Eview;%main form handle
pth = getappdata(evi, 'pth');%directory
flNm = getappdata(evi, 'flNm');%filename

%{'*.eps', 'eps'; '*.emf', 'emf'; '*.jpeg', 'jpeg'}
[imgFlNm, pth] = uiputfile({'*.emf', 'emf'}, 'Save as', [pth, flNm]);%save image dialog
z = find(imgFlNm == '.', 1, 'last');
saveas(evi, [pth, imgFlNm], imgFlNm((z + 1):end))
imgM = getframe(evi);
imwrite(imgM.cdata, [pth, imgFlNm(1:(z - 1)), '.jpeg'], 'jpeg')

% imgM = getframe(handles.LFP_ax);
% imwrite(imgM.cdata, [pth, imgFlNm(1:(z - 1)), '_LFP-CSD', imgFlNm(z:end)], imgFlNm((z + 1):end))
% imgM = getframe(handles.Spk_ax);
% imwrite(imgM.cdata, [pth, imgFlNm(1:(z - 1)), '_MUA', imgFlNm(z:end)], imgFlNm((z + 1):end))

% -print


% --- Executes on button press in SelectChann.
function SelectChann_Callback(hObject, eventdata, handles)
% hObject    handle to SelectChann (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%evi = handles.Eview;%main form handle
chTab = handles.ShowChannels;

if isequal(get(chTab, 'Visible'), 'on')
    set(chTab, 'Visible', 'off')
    set(handles.SelectAllCh, 'Visible', 'off')
    set(handles.InvertSelection, 'Visible', 'off')
    set(handles.SpkShow, 'Visible', 'on')
    PlotAll(handles);
else
    set(chTab, 'Visible', 'on')
    set(handles.SelectAllCh, 'Visible', 'on')
    set(handles.InvertSelection, 'Visible', 'on')
    set(handles.SpkShow, 'Visible', 'off')
end



function MUAScale_Callback(hObject, eventdata, handles)
% hObject    handle to MUAScale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MUAScale as text
%        str2double(get(hObject,'String')) returns contents of MUAScale as a double

PlotAll(handles);


% --- Executes on button press in SelectAllCh.
function SelectAllCh_Callback(hObject, eventdata, handles)
% hObject    handle to SelectAllCh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

chTab = handles.ShowChannels;%table of channels selections
tabData = get(chTab, 'Data');%data of the table
for t = 1:size(tabData, 1) %run over rows of table (channels)
    tabData{t, 2} = true;%select all channels
end
set(chTab, 'Data', tabData)%set new table data



% --- Executes on button press in InvertSelection.
function InvertSelection_Callback(hObject, eventdata, handles)
% hObject    handle to InvertSelection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

chTab = handles.ShowChannels;%table of channels selections
tabData = get(chTab, 'Data');%data of the table
for t = 1:size(tabData, 1) %run over rows of table (channels)
    tabData{t, 2} = ~tabData{t, 2};%invert the selection
end
set(chTab, 'Data', tabData)%set new table data


% --- Executes on button press in SpkShow.
function SpkShow_Callback(hObject, eventdata, handles)
% hObject    handle to SpkShow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

evi = handles.Eview;
lfp_ax = handles.LFP_ax;
spk_ax = handles.Spk_ax;
if get(hObject, 'Value') %checked to show spike axis
    set(lfp_ax, 'Position', [0.05318221447253707, 0.054959785522788206, 0.4272013949433304, 0.8552278820375335])
    set(spk_ax, 'Visible', 'on')
%     chld = get(spk_ax, 'Children');
%     for t = 1:length(chld)
%         set(chld(t), 'Visible', 'on');
%     end
    set(handles.MUAscb, 'Visible', 'on')
    set(handles.Spk_clb, 'Visible', 'on')
    set(handles.MUAmapMin, 'Visible', 'on')
    set(handles.MUAmapMax, 'Visible', 'on')
else %checked to hide spike axis
    set(lfp_ax, 'Position', [0.05318221447253707, 0.054959785522788206, 0.9319965126416739, 0.8552278820375335])
    set(spk_ax, 'Visible', 'off')
%     chld = get(spk_ax, 'Children');
%     for t = 1:length(chld)
%         set(chld(t), 'Visible', 'off');
%     end
    set(handles.MUAscb, 'Visible', 'off')
    set(handles.Spk_clb, 'Visible', 'off')
    set(handles.MUAmapMin, 'Visible', 'off')
    set(handles.MUAmapMax, 'Visible', 'off')
end
PlotAll(handles)


% --- Executes on button press in WinAvr.
function WinAvr_Callback(hObject, eventdata, handles)
% hObject    handle to WinAvr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(hObject, 'Value')
    set(handles.Average, 'Value', 1);
    PlotAll(handles)
end



function MapMinMax_Change(hObject, eventdata, handles)
% hObject    handle to MUAmapMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

PlotAll(handles)


% --------------------------------------------------------------------
function SaveScreenMat_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to SaveScreenMat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%export data of displayed traces to mat-file
lfpTrcs = findobj(get(handles.LFP_ax, 'Children'), 'Type', 'line');
lfpTrcs = findobj(lfpTrcs, 'Color', [0, 0, 0]);
spkTrcs = findobj(get(handles.Spk_ax, 'Children'), 'Type', 'line');
spkTrcs = findobj(spkTrcs, 'Color', [0, 0, 1]);
spkshown = get(handles.SpkShow, 'Value');%is spikes axis visible
if ((length(lfpTrcs) ~= length(spkTrcs)) && spkshown)
    disp('not count')
    return
end
tm_p = get(lfpTrcs(1), 'XData');%lfp time vector
lfp = zeros(length(tm_p), length(lfpTrcs));
if spkshown %spikes visible
    tm_s = get(spkTrcs(1), 'Xdata');%spikes time vector
    spk = zeros(length(tm_s), length(spkTrcs));%spikes density
else
    tm_s = [];
    spk = [];
end
k = 1;
for t = length(lfpTrcs):-1:1 %run over traces
    dCoeffP = get(lfpTrcs(t), 'UserData');
    lfp(:, k) = (get(lfpTrcs(t), 'YData') - dCoeffP{2}) * dCoeffP{1}(3);
    if spkshown %spikes visible
        dCoeffS = get(spkTrcs(t), 'UserData');
        spk(:, k) = (get(spkTrcs(t), 'YData') - dCoeffS{2}) * dCoeffS{1}(3);
    end
    k = k + 1;
end
tm_p = tm_p * dCoeffP{1}(4);%lfp time vector in milliseconds
if spkshown %spikes visible
    tm_s = tm_s * dCoeffS{1}(4);%spikes time vector in milliseconds
end

pth = getappdata(handles.Eview, 'pth');%directory
flNm = getappdata(handles.Eview, 'flNm');%filename
segms = getappdata(handles.Eview, 'segms');
if (get(handles.Average, 'Value') || get(handles.WinAvr, 'Value')) %average segments
    z = str2double(get(handles.CurrSegm, 'String'));%first segment to be averaged
    t = str2double(get(handles.AvrLen, 'String'));%averaging length (number of segments to be averaged)
    sn = max(1, z):min((z + t - 1), size(segms, 1));%selected segments
    %sn = 1:size(segms, 1);%all segments
else %single segment
    sn = str2double(get(handles.CurrSegm, 'String'));%wanted sweeps or segments
    if (sn < 1)
        sn = 1;
    elseif (sn > size(segms, 1))
        sn = size(segms, 1);
    end
end
if (length(sn) == 1)
    bufStr = ['_sn', num2str(sn)];
else
    bufStr = ['_sn', num2str(sn(1)), 'to', num2str(sn(end))];
end
flNm = regexprep(flNm, '_', 'T');
flNm = regexprep(flNm, '-', '_');
lfpName = ['lfp_', flNm, bufStr];%lfp variable name
spkName = ['spk_', flNm, bufStr];%spike varialb name
eval([lfpName, ' = [tm_p'', lfp];']);%lfp values
eval([spkName, ' = [tm_s'', spk];'])%spikes values

if exist('eview_screen.mat', 'file') %storage exist
    save('eview_screen.mat', lfpName, '-append')
else %no storage
    save('eview_screen.mat', lfpName)
end
if eval(['~isempty(', spkName, ')']) %spikes are shown
    save('eview_screen.mat', spkName, '-append')
end


% --- Executes when entered data in editable cell(s) in ShowChannels.
function ShowChannels_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to ShowChannels (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)

%if grpSel edited (change group of cells)

if (eventdata.Indices(2) == 3) %group selection
    tabData = get(handles.ShowChannels, 'Data');
    hd = getappdata(handles.Eview, 'hd');
    z = eventdata.Indices(1);%selected row
    for t = max(1, z - 4):min(z + 4, hd.nADCNumChannels) %run over selected channels
        tabData{t, 2} = ~tabData{t, 2};%invert the selection of the channel (%true;%select group of channels)
    end
    set(handles.ShowChannels, 'Data', tabData)%set table data
end
