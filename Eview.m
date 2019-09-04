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

% Last Modified by GUIDE v2.5 14-Jun-2019 17:55:48

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
% handles.LFPscb = annotation('textarrow', 'Position', [0.025, 0.15 0, 0.1], 'LineWidth', 5, 'Units', 'normalized', 'TextRotation', 90, 'HeadStyle', 'none',...
%     'String', '0 \muV', 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');%LFP-scale bar
% handles.MUAscb = annotation('textarrow', 'Position', [0.02, 0.35 0, 0.1], 'LineWidth', 5, 'Units', 'normalized', 'TextRotation', 90, 'HeadStyle', 'none',...
%     'String', '0 spikes', 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');%MUA-scale bar
setappdata(evi, 'dataCurs', []);

% handles.CSD_clb = colorbar('peer', handles.LFP_ax, 'Location', 'manual', 'Units', 'normalized', ... %CSD-colorbar
%     'Position', [0.01, 0.71, 0.01, 0.17], 'YAxisLocation', 'left', 'FontSize', 6, 'YTickLabel', []);%CSD-colorbar
% handles.Spk_clb = colorbar('peer', handles.LFP_ax, 'Location', 'manual', 'Units', 'normalized', ... %MUA-colorbar
%     'Position', [0.015, 0.12, 0.01, 0.17], 'YAxisLocation', 'left', 'FontSize', 6, 'YTickLabel', []);%MUA-colorbar

set(evi, 'WindowButtonDownFcn', {@Eview_WindowButtonDownFcn, handles})%reaction on button down
set(evi, 'WindowButtonMotionFcn', {@Eview_WindowButtonMotionFcn, handles})%reaction on button motion
set(evi, 'WindowButtonUpFcn', {@Eview_WindowButtonUpFcn, handles})%reaction on button up

h=uicontextmenu;%contex
hrmb=uimenu(h, 'Label', 'individual', 'Callback', {@SeparateChannel, handles});
set(handles.ChannSelector, 'UIContextMenu', h);%associate context menu with channel-seletor table

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

if ~isempty(getappdata(evi, 'brst')) && get(handles.SaveBursts, 'UserData') %changes was done
    ynBtn = questdlg('Save resultes?', 'Save on exit', 'Yes', 'No', 'No');%ask what to do
    if isequal(ynBtn, 'Yes') %save requested
        SaveBursts_ClickedCallback([], [], handles);%save bursts
    end
end
        
[flNm, pth] = uigetfile({'*.mat', 'mat'}, 'Select file', [pth, flNm, '.mat']);%open dialog
if (isequal(flNm, 0) && isequal(pth, 0)) %no file choosed
    return
end

z = find(flNm == '.', 1, 'last');
if ~isequal(flNm(z:end), '.mat')
    %set(nts, 'String', 'unknown file type')
    return
end
flNm = flNm(1:(z - 1));%delete expansion

varInMatFile = who(matfile([pth, flNm, '.mat']));
if ismember('zp', varInMatFile)
    load([pth, flNm, '.mat'], 'lfp', 'spks', 'hd', 'zp', 'lfpVar')
    zavp = zp;
else
    load([pth, flNm, '.mat'], 'lfp', 'spks', 'hd', 'zavp', 'lfpVar')
end
set(evi, 'Name', ['Eview - ', flNm])
if ~isfield(zavp, 'dwnSmplFrq') %field exist
    zavp.dwnSmplFrq = 1e3;%discretization frequency for signals to be treated (Hz)
end

nlxVer = 0;%old format of Cheetah log-file
if exist([pth, flNm], 'dir') %neuralynx directory
    dirCnt = dir([pth, flNm]);%directory content
    for t = 3:length(dirCnt) %run over files in NLX directory
        if ((length(dirCnt(t).name) > 3) && isequal(dirCnt(t).name(end - 3:end), '.nde'))
            nlxVer = 1;%new format of Cheetah log-file
            break;%out of (for t)
        end
    end
end

%lfp = double(lfp);
prg = str2double(get(handles.SpkPrg, 'String'));
spksOrig = spks;%copy of spikes
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

%= load events from xls =%
z = find(flNm == '_', 1, 'first');
if exist([pth, flNm(1:z), 'events.xlsx'], 'file')
    [~, ~, eventsTags] = xlsread([pth, flNm(1:z), 'events.xlsx']);%load manual protocol
    for t = 2:length(eventsTags) %run over events
        eventsTags{t} = ((eventsTags{t} * 3600 * 24) - hd.recTime(1)) * 1e3;%convert to ms from current record begin
    end
else
    eventsTags = {};
end

%= detected bursts (any segments of recordation) =%
brst = [];%clear previous bursts
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
    set(handles.ShowBrst, 'Enable', 'on')
end
set(handles.SaveBursts, 'UserData', 0)%set flag of bursts saved

%= fill table of channel selector =%
chTabData = get(handles.ChannSelector, 'Data');%channel selection table
if isempty(horzcat(chTabData{:, 1})) %first initiation
    chTabData = cell(hd.nADCNumChannels, 2);%create channels-table data
    
    %creat channel groups with default settings 
    chanSettingGrp(1:hd.nADCNumChannels) = struct('grpName', [], ... name of channels group
                                                  'lowCut', false, ... do not lowcut filtering
                                                  'lowCutFrq', 0.1, ... lowcut frequency
                                                  'highCut', false, ... do not highcut filtering
                                                  'highCutFrq', 100, ... highcut frequency
                                                  'rawSignal', false, ... raw signal
                                                  'comRef', false, ... common reference
                                                  'dcShift', false, ... compensate DC
                                                  'pltCSDmap', false, ... plot CSD map
                                                  'pltSpkMap', false, ... plot spikes frequency map
                                                  'pltLFP', true, ... plot LFP traces
                                                  'pltMUA', true, ... plot MUA lines
                                                  'pltMUAfreq', false, ... plot trace of spikes densities
                                                  'scale', 1 ... verticale scale (mV, amplitude)
                                                  );
    for t = 1:hd.nADCNumChannels %run over channels
        chTabData{t, 3} = false;%channel deselected
        chanSettingGrp(t).grpName = t;%initially each channel is in its own group
        chTabData{t, 4} = chanSettingGrp(t).grpName;%copy group name
    end
else %get previous settings
    chanSettingGrp = getappdata(evi, 'chanSettingGrp');%load channels group setting
    if (size(chTabData, 1) >= hd.nADCNumChannels) %current file contain less channels
        chTabData((hd.nADCNumChannels + 1):end, :) = [];%delete bottom part of channels list
        chanSettingGrp((hd.nADCNumChannels + 1):end) = [];%delete bottom part of settings list
    else %current file contain more channels (creat setting for additional channels)
        z = size(chTabData, 1);%number of last previous channel
        chTabData((z + 1):hd.nADCNumChannels, :) = cell((hd.nADCNumChannels - z), 4);%add cells
        chanSettingGrp((z + 1):hd.nADCNumChannels) = struct('grpName', [], ... name of channels group
                                                  'lowCut', false, ... do not lowcut filtering
                                                  'lowCutFrq', 0.1, ... lowcut frequency
                                                  'highCut', false, ... do not highcut filtering
                                                  'highCutFrq', 100, ... highcut frequency
                                                  'rawSignal', false, ... raw signal
                                                  'comRef', false, ... common reference
                                                  'dcShift', false, ... compensate DC
                                                  'pltCSDmap', false, ... plot CSD map
                                                  'pltSpkMap', false, ... plot spikes frequency map
                                                  'pltLFP', true, ... plot LFP traces
                                                  'pltMUA', true, ... plot MUA lines
                                                  'pltMUAfreq', false, ... plot trace of spikes densities
                                                  'scale', 1 ... verticale scale (mV, amplitude)
                                                  );%add structures
        m = max(horzcat(chanSettingGrp(1:z).grpName)) + 1;%new groups for each new channel
        for t = (z + 1):hd.nADCNumChannels %run over channels
            chTabData{t, 3} = false;%channel deselected
            chanSettingGrp(t).grpName = m;%each new channel is in its own group
            chTabData{t, 4} = chanSettingGrp(t).grpName;%copy group name
            m = m + 1;%new group for each new channel
        end
    end
end

%rewrite names of channels
if isfield(hd, 'recChNames') %channel names exist
    for t = 1:hd.nADCNumChannels %run over channels
        chTabData{t, 1} = hd.recChNames{t, 1};%name of channel
        chTabData{t, 2} = '';%altername of channel
    end
    if (size(hd.recChNames, 2) > 1) %altername was signed
        for t = 1:hd.nADCNumChannels %run over channels
            chTabData{t, 2} = hd.recChNames{t, 2};%name of channel
        end
    end
else %channel names not exist
    hd.recChNames = cell(hd.nADCNumChannels, 1);
    for t = 1:hd.nADCNumChannels %run over channels
        chTabData{t, 1} = num2str(t);%name of channel
    end
end
if ~any(vertcat(chTabData{:, 3})) %previous selection
    chTabData{1, 3} = true;%select only first channel
end
set(handles.ChannSelector, 'Data', chTabData)%set table data

%delete "NaN" and "Inf" from lfp
for ch = 1:size(lfp, 2) %run over channels
    p1 = find(~isfinite(lfp(:, ch)), 1, 'first');
    while ~isempty(p1)
        p2 = p1 + find(isfinite(lfp(p1:end, ch)), 1, 'first') - 2;
        lfp(p1:p2, ch) = interp1([p1; p2], [nanmedian(lfp(p1 + (-100:-1), ch)), nanmedian(lfp(p2 + (1:100), ch))], (p1:p2)');
        p1 = find(~isfinite(lfp(:, ch)), 1, 'first');
    end
end

%= save application data =%
setappdata(evi, 'pth', pth) %directory
setappdata(evi, 'flNm', flNm) %filename
setappdata(evi, 'lfp', lfp) %LFP
setappdata(evi, 'spks', spks) %separated spikes
setappdata(evi, 'spksOrig', spksOrig) %original spikes
setappdata(evi, 'lfpVar', lfpVar) %LFP variance
setappdata(evi, 'hd', hd) %header
setappdata(evi, 'zavp', zavp) %additional parameters
setappdata(evi, 'eventsTags', eventsTags) %manually entered events
setappdata(evi, 'brst', brst) %bursts
setappdata(evi, 'selectedBrst', []) %selected burst
setappdata(evi, 'nlxVer', nlxVer) %Neuralynx cheetah version
setappdata(evi, 'chanSettingGrp', chanSettingGrp);%save channels group setting

CutSpnt_Callback([], [], handles) %creat segment and plot


function BefStim_Callback(hObject, eventdata, handles)
% hObject    handle to BefStim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

PlotAll(handles);%plot with new parameters


function AftStim_Callback(hObject, eventdata, handles)
% hObject    handle to AftStim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

PlotAll(handles);%plot with new parameters



% --- Executes on button press in NSegm.
function PNSegm_Callback(hObject, eventdata, handles)
% hObject    handle to NSegm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

evi = handles.Eview;%handle of the application window
crSeg = handles.CurrSegm;
segms = getappdata(evi, 'segms');%number of signal segments
cs = str2double(get(crSeg, 'String'));%current segment
if isequal(hObject, handles.NSegm)
    cs = cs + str2double(get(handles.Step, 'String'));
elseif isequal(hObject, handles.PSegm)
    cs = cs - str2double(get(handles.Step, 'String'));
end
if (cs < 1)
    cs = 1;
elseif (cs > size(segms, 1))
    cs = size(segms, 1);
end
set(crSeg, 'String', num2str(cs))
PlotAll(handles);%plot with new parameters



function CurrSegm_Callback(hObject, eventdata, handles)
% hObject    handle to CurrSegm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

PlotAll(handles);%plot with new parameters



% --- Executes on button press in CutSpnt.
function CutSpnt_Callback(hObject, eventdata, handles)
% hObject    handle to CutSpnt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
evi = handles.Eview;%handle of the application window
hd = getappdata(evi, 'hd');
lfp = getappdata(evi, 'lfp');
zavp = getappdata(evi, 'zavp');

%= creat segments of the recordation =%
%segms(:, 1) - position of synchro-event (ms)
%segms(:, 2) - number of channel where syncro-event was detected
%segms(:, 3) - number of sweeps where syncro-event was detected
%segms(:, 4) - number of stimulus (within sweep, in matrix zavp.realStim) or number of trough in brst matrix
if (get(handles.CutSpnt, 'Value') || (numel(zavp.realStim(1).r) <= 0)) %spontwise cutting choosed
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
else %cut with respect to synchronization
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
    tW = min(1e5, round(size(lfp, 1) ./ size(segms, 1)));%length of view window (ms)
end

%= creat time-tags =%
timeStr = repmat(' ', size(segms, 1), 8);%"Moscow" time of the segment begin (part of day)
for t = 1:size(segms, 1) %run over segments
    timeStr(t, :) = MoscowTime(hd.recTime(1) + (segms(t, 1) * 1e-3));%"Moscow" time of the segment begin (part of day)
end

%= save application data =%
setappdata(evi, 'segms', segms) %segments of recordation
setappdata(evi, 'timeStr', timeStr) %automatically created time tags

set(handles.CurrSegm, 'String', '1')
set(handles.SegmTxt, 'String', num2str(size(segms, 1)))
if (isnan(str2double(get(handles.BefStim, 'String'))) || isnan(str2double(get(handles.AftStim, 'String'))))
    set(handles.BefStim, 'String', '0')
    set(handles.AftStim, 'String', num2str(tW))%length of view window
end

PlotAll(handles);%plot with new parameters



function PlotAll(handles)
evi = handles.Eview;%handle of the application window
lfp_ax = handles.LFP_ax;%LFP axis handle

dataCurs = getappdata(handles.Eview, 'dataCurs');%data cursors
for t = 1:length(dataCurs)
    delete(dataCurs(t).t);%delete data cursor
    delete(dataCurs(t).l);%delete line
end
setappdata(handles.Eview, 'dataCurs', []);%data cursors

%= load application data =%
% pth = getappdata(evi, 'pth');%directory
flNm = getappdata(evi, 'flNm');%filename
lfp = getappdata(evi, 'lfp');%LFP
spks = getappdata(evi, 'spks');%separated spikes
% spksOrig = getappdata(evi, 'spksOrig');%original spikes
lfpVar = getappdata(evi, 'lfpVar');%LFP variance
hd = getappdata(evi, 'hd');%header
zavp = getappdata(evi, 'zavp');%additional parameters
eventsTags = getappdata(evi, 'eventsTags');%manually entered events
brst = getappdata(evi, 'brst');%bursts
% selectedBrst = getappdata(evi, 'selectedBrst');%selected burst
nlxVer = getappdata(evi, 'nlxVer');%Neuralynx cheetah version
chanSettingGrp = getappdata(evi, 'chanSettingGrp');%save channels group setting
chTabData = get(handles.ChannSelector, 'Data');%channels selector table content
segms = getappdata(evi, 'segms');%segments of recordation
timeStr = getappdata(evi, 'timeStr');%automatically created time tags

if isfield(zavp, 'chOrd')
    lfp = lfp(:, zavp.chOrd);
    spks = spks(zavp.chOrd, :);
    lfpVar = lfpVar(zavp.chOrd, 1);%LFP variance (with respect to daq-files)
end

segmEdge = [-str2double(get(handles.BefStim, 'String')), ... left shifts from synchro-point (ms)
             str2double(get(handles.AftStim, 'String'))];%right shifts from synchro-point (ms)

if (get(handles.Average, 'Value') || get(handles.WinAvr, 'Value')) %average segments
    z = str2double(get(handles.CurrSegm, 'String'));%first segment to be averaged
    t = str2double(get(handles.AvrLen, 'String'));%averaging length (number of segments to be averaged)
    sn = max(1, z):min((z + t - 1), size(segms, 1));%selected segments
    %sn = 1:size(segms, 1);%all segments
else %single segment
    sn = str2double(get(handles.CurrSegm, 'String'));%wanted sweeps or segments
    if (sn < 1)
        sn = 1;%correct number
    elseif (sn > size(segms, 1))
        sn = size(segms, 1);%correct number
    end
end
set(handles.CurrSegm, 'String', num2str(sn(1)));%show number of of current segment
set(handles.Average, 'Value', 0);%next plot single (change flag)

%= get numbers (names) of unique group =%


cla(lfp_ax)%clear LFP-axis
setappdata(evi, 'dFrmH', []);%drag frame handle
setappdata(evi, 'bgnXY', []);%begin point of rectangle draw
setappdata(evi, 'endXY', []);%end point of rectangle draw
chnlStr = cell(0);%vertical axis labels (channels)

%= run over channel groups and plot traces individual settings =%
unqChnlGrp = unique(horzcat(chTabData{:, 4}));%names of unique groups
initLvl = 0;%initial drawing level
for cg = unqChnlGrp %run over unique groups
    jj = horzcat(chanSettingGrp(:).grpName);%all groups
    cgch = find(jj == cg);%current group channels
    rCh = find(horzcat(chTabData{:, 3}));%all selected channels
    rCh = intersect(rCh, cgch);%number of visible (selected) channels of current group rCh array
    
    if ~isempty(rCh) %any channels was selected
    rawData = chanSettingGrp(rCh(1)).rawSignal;%load raw data if true
    
    %= timevector for raw and resampled signals =%
    if rawData %raw data requested
        prm.Fs = 1 / zavp.siS;%sampling frequency for raw data
    else %resampled data requested
        prm.Fs = zavp.dwnSmplFrq;%sampling frequency for rare data
    end
    tm = segmEdge(1):(1e3 / prm.Fs):segmEdge(2);%time matrix for raw or downsampled data (ms)

    %= common reference (step 1) =%
    if chanSettingGrp(rCh(1)).comRef %common reference subtraction
        fLFP = lfp(:, cgch);%no filters
        %= filtration =%
        if (chanSettingGrp(rCh(1)).lowCut && (chanSettingGrp(rCh(1)).lowCutFrq < (prm.Fs / 2))) %lowcut filtration
            if (chanSettingGrp(rCh(1)).lowCutFrq < 2) %small cutting frequency
                fLFP = ZavRCfilt(lfp(:, cgch), chanSettingGrp(rCh(1)).lowCutFrq, prm.Fs, 'high');%RC-filter, lowcut
            else %quite large cutting frequency
                fLFP = ZavFilter(double(lfp(:, cgch)), prm.Fs, 'high', chanSettingGrp(rCh(1)).lowCutFrq, 2);%lowcut
            end
        end
        if (chanSettingGrp(rCh(1)).highCut && (chanSettingGrp(rCh(1)).highCutFrq < (prm.Fs / 2))) %highcut filtration
            if (chanSettingGrp(rCh(1)).highCutFrq < 2) %small cutting frequency
                fLFP = ZavRCfilt(lfp(:, cgch), chanSettingGrp(rCh(1)).highCutFrq, prm.Fs, 'low');%RC-filter, lowcut
            else
                fLFP = ZavFilter(double(lfp(:, cgch)), prm.Fs, 'low', chanSettingGrp(rCh(1)).highCutFrq, 2);%highcut
            end
        end
    
        whtLFP = ZavSynchLFP(zavp, hd, segms(sn, :), segmEdge, fLFP, 1:size(fLFP, 2), rawData, nlxVer);%LFP phased with respect to stimulus moments
    else %without common reference subtraction
        fLFP = lfp(:, rCh);%no filters
        %= filtration =%
        if (chanSettingGrp(rCh(1)).lowCut && (chanSettingGrp(rCh(1)).lowCutFrq < (prm.Fs / 2))) %lowcut filtration
            if (chanSettingGrp(rCh(1)).lowCutFrq < 2) %small cutting frequency
                fLFP = ZavRCfilt(lfp(:, rCh), chanSettingGrp(rCh(1)).lowCutFrq, prm.Fs, 'high');%RC-filter, lowcut
            else %quite large cutting frequency
                fLFP = ZavFilter(double(lfp(:, rCh)), prm.Fs, 'high', chanSettingGrp(rCh(1)).lowCutFrq, 2);%lowcut
            end
        end
        if (chanSettingGrp(rCh(1)).highCut && (chanSettingGrp(rCh(1)).highCutFrq < (prm.Fs / 2))) %highcut filtration
            if (chanSettingGrp(rCh(1)).highCutFrq < 2) %small cutting frequency
                fLFP = ZavRCfilt(lfp(:, rCh), chanSettingGrp(rCh(1)).highCutFrq, prm.Fs, 'low');%RC-filter, lowcut
            else
                fLFP = ZavFilter(double(lfp(:, rCh)), prm.Fs, 'low', chanSettingGrp(rCh(1)).highCutFrq, 2);%highcut
            end
        end
        
        whtLFP = ZavSynchLFP(zavp, hd, segms(sn, :), segmEdge, fLFP, 1:size(fLFP, 2), rawData, nlxVer);%lfp phased with respect to stimulus moments
    end

    %= DC-shift compensation =%
    if chanSettingGrp(rCh(1)).dcShift %DC compensation
        %ii = min(5, size(whtLFP, 1)):min(1e3, size(whtLFP, 1));
        ii = (tm > -1e1) & (tm < 1e1);%time near synchro-event
        ajj = median(whtLFP(ii, :, :), 1);%base-level
        whtLFP = whtLFP - repmat(ajj, [size(whtLFP, 1), 1, 1]);%subtract DC component
    end

    %= common reference (step 2) =%
    if chanSettingGrp(rCh(1)).comRef %common reference subtraction
        %ii = min(5, size(whtLFP, 1)):min(1e3, size(whtLFP, 1));
        ii = (tm > -1e1) & (tm < 1e1);%time near synchro-event
        ajj = median(whtLFP(ii, :, :), 1);%base-level
        wLFPcomm = mean(whtLFP - repmat(ajj, [size(whtLFP, 1), 1, 1]), 2);%common reference with adduction to zero-level of each channel in a group
        [~, irch] = intersect(cgch, rCh);%number of selected (visible) channel of the group
        whtLFP = whtLFP(:, irch, :) - repmat(wLFPcomm, [1, length(rCh), 1]);%subtract common reference
    end
    whtLFP = mean(whtLFP, 3);%average by segments (after transformations)
    

    %= time axis parameters =%
    if (diff(segmEdge) > 300e3) %minutes
        timeScalCoeff = 60e3;%time coefficient
        timeUnit = 'min';%units
    elseif ((diff(segmEdge) > 5e3) && (diff(segmEdge) <= 300e3)) %seconds
        timeScalCoeff = 1e3;%time coefficient
        timeUnit = 's';%units
    else %milliseconds
        timeScalCoeff = 1;%time coefficient
        timeUnit = 'ms';%units
    end
    
    %= plot with step (drawing acceleration) =%
    rarStep = 1 + floor(size(whtLFP, 1) / 10e3);%plot with step
    if (rarStep > 1)
        ajj = whtLFP;%copy of LFP segment
        whtLFP = zeros(length(resample(ajj(:, 1), 1, rarStep)), size(ajj, 2));
        for t = 1:size(whtLFP, 2) %run over channels
            whtLFP(:, t) = resample(ajj(:, t), 1, rarStep);%downsampled LFP
        end
        tm = interp1(1:length(tm), tm, linspace(1, length(tm), size(whtLFP, 1)));%downsampled time
    end

    %= plot CSD =%
    if (chanSettingGrp(rCh(1)).pltCSDmap && (size(whtLFP, 2) > 2)) %plot current source density (CSD)
        if (size(whtLFP, 1) > 1e3) %long segment
            ajj = ZavRCfilt(whtLFP, 0.5, 1e3, 'high');%remove DC
        else %short segmen
            ajj = whtLFP - repmat(mean(whtLFP(1:min(10, size(whtLFP, 1)), :), 1), length(tm), 1);%remove DC
        end
    %     ajj = ajj ./ repmat(lfpVar(rCh, 1)' ./ min(lfpVar(rCh, 1)), length(tm), 1);%amplitude correction

        lfpCSD = -diff(ajj, 2, 2);%CSD
        k = 3;%round(diff(segmEdge) / 100) + 1;%interval of smoothing
        for t = 1:size(lfpCSD, 2) %run over channels
            lfpCSD(:, t) = smooth(lfpCSD(:, t), k);%smoothing along time
        end
        %interpolation (along channels)
        ajj = lfpCSD;%copy of original CSD map
        lfpCSD = zeros(size(ajj, 1), length(1:0.2:size(lfpCSD, 2)));
        for t = 1:size(lfpCSD, 1) %run over time
            lfpCSD(t, :) = interp1(1:(size(whtLFP, 2) - 2), ajj(t, :), 1:0.2:(size(whtLFP, 2) - 2));%interpolation along channels
        end

        subplot(lfp_ax);
        imagesc(tm / timeScalCoeff, initLvl + (2:(length(rCh) - 1)), lfpCSD');
        caxis(lfp_ax, [str2double(get(handles.CSDcMapMin, 'String')) / 1e3, str2double(get(handles.CSDcMapMax, 'String')) / 1e3])%colorbar range
    end

    %= plot spikes density =%
    if (chanSettingGrp(rCh(1)).pltSpkMap || chanSettingGrp(rCh(1)).pltMUAfreq) %plot map of spikes frequency
        binStep = [str2double(get(handles.SpkBin, 'String')), str2double(get(handles.HistPrecision, 'String'))];%[bin size (ms), precision (must be integer)]
        [~, basbins] = ZavOverHist(spks(rCh(1), 1).tStamp, binStep, segmEdge);%fast overlaped histogram (histc)
        spkDens = zeros(length(basbins), numel(rCh), numel(sn));%memory preallocation for density of spikes in time (histogram)

        for z = 1:numel(sn) %run over requested segments
            jj = segms(sn(z), 1) + segmEdge;%time segment boundary (ms from begin of file)
            sw = segms(sn(z), 3);%sweep
            for ch = 1:numel(rCh) %run over channels
                spkDens(:, ch, z) = ZavOverHist(spks(rCh(ch), sw).tStamp, binStep, jj);%fast overlaped histogram (histc)
            end
        end
        spkDens = mean(spkDens, 3) / binStep(1);%mean spikes density (spikes / ms)

        p1 = str2double(get(handles.MUAScale, 'String'));
        if ~isempty(p1)
            if ~isfinite(p1)
                p1 = max(max(spkDens));
                set(handles.MUAScale, 'String', num2str(p1))
            end
        end

        dCoeff(1) = -1; dCoeff(2) = p1;%characteristic sizes
        dCoeff(3) = dCoeff(1) * dCoeff(2);%scaling coefficient

        clr = 'b';%color of spikes density line
        if chanSettingGrp(rCh(1)).pltSpkMap %plot MUA density map
        %     %interpolation (smooth along channels)
        %     spkDensSC = zeros(size(spkDens, 1), length(1:0.25:length(rCh)));
        %     for t = 1:size(spkDens, 1) %run over time
        %         spkDensSC(t, :) = interp1(1:length(rCh), spkDens(t, :), 1:0.25:length(rCh));%interpolation along channels
        %     end
            spkDensSC = spkDens;
            %smooth along time
            for t = 1:size(spkDens, 2) %rum over channels
                spkDensSC(:, t) = smooth(spkDens(:, t), 3);%smooth along time
            end

            imagesc(basbins / timeScalCoeff, initLvl + (1:length(rCh)), spkDensSC');
            shading flat;
            caxis(lfp_ax, [str2double(get(handles.SpkCMapMin, 'String')), str2double(get(handles.SpkCMapMax, 'String'))])%colorbar range
            clr = 'w';%color of spikes density line
        end %if ~chanSettingGrp(rCh(1)).pltCSDmap %no CSD plot choosed

        %= plot spikes densities
        if chanSettingGrp(rCh(1)).pltMUAfreq %plot MUA frequency traces
            linHandl = plot(lfp_ax, basbins / timeScalCoeff, (spkDens / dCoeff(3)) + repmat(initLvl + (1:length(rCh)), length(basbins), 1), clr, 'LineWidth', 1);
            set(linHandl, 'ButtonDownFcn', {@Line_ButtonDownFcn, handles})%reaction on mouse-click
            for t = 1:length(linHandl) %run over lines handles
                set(linHandl(t), 'UserData', {dCoeff, t, timeScalCoeff});%save amplitude coefficient and offset
            end

    %         [~, p1] = ds2nfu([0, 0], [0, abs(1 / dCoeff(1))]);%diff(p1) = length of scale bar in normalized units
    %         set(handles.MUAscb, 'Position', [0.525, 0.15 0, diff(p1)], 'LineWidth', 5, 'Units', 'normalized', 'TextRotation', 90, 'HeadStyle', 'none',...
    %             'String', [num2str(round(dCoeff(2) * 10) / 10), 'spikes/', num2str(binStep(1)), 'ms'], 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');
        end
    end
    
    %= plot LFP =%
    if chanSettingGrp(rCh(1)).pltLFP %plot LFP traces
        %= LFP scaling =%
        dCoeff(1) = -1.0; dCoeff(2) = chanSettingGrp(rCh(1)).scale;%characteristic sizes (amplitude coefficient)
        dCoeff(3) = dCoeff(1) * dCoeff(2) * 1e3;%scaling coefficient
        whtLFP = whtLFP / dCoeff(3);%scaled signals
    
        linHandl = plot(lfp_ax, tm / timeScalCoeff, whtLFP + repmat(initLvl + (1:length(rCh)), length(tm), 1), 'k', 'LineWidth', 1);
        set(linHandl, 'ButtonDownFcn', {@Line_ButtonDownFcn, handles})%reaction on mouse-click
        for t = 1:length(linHandl) %run over lines handles
            set(linHandl(t), 'UserData', {dCoeff, t, timeScalCoeff});%save amplitude coefficient and offset
        end

        %= plot bursts =%
        stt = get(handles.ShowBrst, 'State');%plot bursts?
        if (~isempty(brst) && (length(sn) == 1) && isequal(stt, 'on')) %single segment
            jj = vertcat(brst(:).t);
            p1 = find((jj(:, 1) < (segms(sn, 1) + segmEdge(2))) & (jj(:, 2) > (segms(sn, 1) + segmEdge(1))));%find visible bursts
            for t = p1' %run over visible bursts
                [~, ii] = intersect(rCh, brst(t).ch);%common channels
                if ~isempty(ii) %some channels contain the burst
                    jj = ((tm >= (brst(t).t(1) - segms(sn, 1))) & (tm <= (brst(t).t(2) - segms(sn, 1))));
                    plot(lfp_ax, tm(jj) / timeScalCoeff, whtLFP(jj, ii) + repmat(initLvl + (rCh(ii) - (rCh(1) - 1)), sum(jj), 1), 'm', 'LineWidth', 1)
                end
            end
        end
    end

    %= plot spikes =%
    if ((length(sn) == 1) && chanSettingGrp(rCh(1)).pltMUA) %single sweep (segment)
        jj = segms(sn, 1) + segmEdge;%absolute time of the segment boundary ( ms, e.i. from begin of file)
        sw = segms(sn, 3);%sweep
        for ch = 1:length(rCh) %run over wanted channels
            spkDens = spks(rCh(ch), sw).tStamp((spks(rCh(ch), sw).tStamp > jj(1)) & (spks(rCh(ch), sw).tStamp < jj(end))) - segms(sn, 1);%unites in the segment
            spkCdens = NaN(3 * length(spkDens), 2);%continuous trace of spikes
            spkCdens(1:3:end, 1) = spkDens;%time point of lower level
            spkCdens(2:3:end, 1) = spkDens;%time point of upper level
            spkCdens(1:3:end, 2) = initLvl + ch - 0.4;%lower level
            spkCdens(2:3:end, 2) = initLvl + ch - 0.2;%upper level
            %plot(lfp_ax, [spkDens'; spkDens'] / timeScalCoeff, t - repmat([0.4; 0.1], 1, length(spkDens)), 'r', 'LineWidth', 1.5);
            plot(lfp_ax, spkCdens(:, 1) / timeScalCoeff, spkCdens(:, 2), 'r', 'LineWidth', 1.5)
        end
    end

    % p1 = str2double(get(handles.LFPScale, 'String'));
    % if ~isempty(p1)
    %     if ~isfinite(p1)
    %         p1 = max(max(abs(whtLFP)));
    %         set(handles.LFPScale, 'String', num2str(p1))
    %     end
    % end

    %= vertical axis labels =%
    if (size(hd.recChNames, 2) > 1) %alternames was signed
        chnlStr = [chnlStr; hd.recChNames(rCh, 2)];%alter name of channels
    else
        chnlStr = [chnlStr; num2cell(rCh(:))];%numbered channels
    end

    % set(get(lfp_ax, 'YLabel'), 'String', 'channels', 'Units', 'normalized', 'Position', [-0.06, 0.6, 1])
    % [~, p1] = ds2nfu([0, 0], [0, abs(1 / dCoeff(1))]);%diff(p1) = length of scale bar in normalized units
    % set(handles.LFPscb, 'Position', [0.025, 0.02, 0, diff(p1)], 'LineWidth', 5, 'Units', 'normalized', 'TextRotation', 90, 'HeadStyle', 'none',...
    %     'String', [num2str(round(dCoeff(2))), ' \muV'], 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', 'FontSize', 9);
    
    initLvl = initLvl + length(rCh);%drawing level for next group
    end %end of (if ~isempty(rCh) %any channels was selected)
end %end of (for cg = unqChnlGrp %run over unique groups)

%= complete axis settings =%
rCh = find(horzcat(chTabData{:, 3}));%all selected channels
y_lmt = [-0.05, length(rCh) + 1.05];%vertical limits
set(lfp_ax, 'YDir', 'reverse', 'XLim', segmEdge / timeScalCoeff, 'YLim', y_lmt)
set(get(lfp_ax, 'XLabel'), 'String', ['time, ', timeUnit])
set(lfp_ax, 'YTick', 1:length(chnlStr), 'YTickLabel', chnlStr)
plot(lfp_ax, [0, 0], y_lmt, 'c', 'LineWidth', 2)%nul-time line (synchro-event)

if (length(sn) == 1) %single segment is shown
    set(evi, 'Name', ['Eview - ', flNm, ' (', timeStr(sn, :), ')'])
    text(0, y_lmt(1), timeStr(sn, :))
    for t = 2:length(eventsTags) %run over events
        if ((eventsTags{t, 1} > (segms(sn, 1) + segmEdge(1))) && (eventsTags{t, 1} < (segms(sn, 1) + segmEdge(2))))
            plot(lfp_ax, (eventsTags{t, 1} - segms(sn, 1)) * [1, 1] / timeScalCoeff, y_lmt, 'g', 'LineWidth', 2) %event-line
            text((eventsTags{t, 1} - segms(sn, 1)) / timeScalCoeff, y_lmt(1), eventsTags{t, 2})
        end
    end
else
    set(evi, 'Name', ['Eview - ', flNm, ' (', timeStr(sn(1), :), ' - ', timeStr(sn(end), :), ')'])
end



% --- Executes on button press in Average.
function Average_Callback(hObject, eventdata, handles)
% hObject    handle to Average (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(hObject, 'Value', 1);
PlotAll(handles);%plot with new parameters



function LFPScale_Callback(hObject, eventdata, handles)
% hObject    handle to LFPScale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

PlotAll(handles);%plot with new parameters



function SpkPrg_Callback(hObject, eventdata, handles)
% hObject    handle to SpkPrg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

evi = handles.Eview;%handle of the application window
prg = str2double(get(handles.SpkPrg, 'String'));%threshold multiplier
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
PlotAll(handles);%plot with new parameters

function timeStr = MoscowTime(timeNum)
%convert time as number to string in format hh:mm
%INPUTS
%timeNum - time as number (seconds from day beginning)
%OUTPUTS
%timeStr - time as string in format hh:mm:ss

tm = zeros(1, 2);
tm(1) = floor(timeNum / 3600);%hours
tm(2) = floor((timeNum - (tm(1) * 3600)) / 60);%minutes
tm(3) = round(timeNum - (tm(1) * 3600) - tm(2) * 60);%seconds

timeStr = num2str(tm(1));
if (tm(1) < 10)
    timeStr = ['0', timeStr];%add zeros
end
timeStr = [timeStr, ':'];
if (tm(2) < 10)
    timeStr = [timeStr, '0'];%add zeros
end
timeStr = [timeStr, num2str(tm(2))];%complete string
timeStr = [timeStr, ':'];
if (tm(3) < 10)
    timeStr = [timeStr, '0'];%add zeros
end
timeStr = [timeStr, num2str(tm(3))];%complete string


function SpkBin_Callback(hObject, eventdata, handles)
% hObject    handle to SpkBin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

PlotAll(handles);%plot with new parameters


% --- Executes on mouse press over axes background.
function Line_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to LFP_ax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

evi = handles.Eview;%handle of the application window
dataCurs = getappdata(evi, 'dataCurs');%data cursors
segms = getappdata(evi, 'segms');%segments of recordation
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
dCoeff = get(hObject, 'UserData');%{amplitude coefficient, number of trace, time coefficient}
timeScalCoeff = dCoeff{3};%time scale coefficient
xy(2) = dCoeff{1}(3) * (xy(2) - dCoeff{2});%signal recovery

cs = str2double(get(handles.CurrSegm, 'String'));%current segment
[~, m] = min(abs(segms(:, 1) - (segms(cs, 1) + xy(1) * timeScalCoeff)));%number of the nearest timestamp

dataCurs(nCrs).t = annotation('textbox', 'Units', 'normalized', 'Position', [p1, 1 - p2 - 0.02, 0.05, 0.05], ...
    'String', {['X:', num2str(round(xy(1) * 10) / 10)], ['Y:', num2str(round(xy(2) * 10) / 10)], timeStr(m, :)}, ...
    'BackgroundColor', [1, 1, 1], 'EdgeColor', [0, 1, 0], 'UserData', {xy, hObject}, 'FontSize', 8);
dataCurs(nCrs).l = annotation('line', 'Units', 'normalized', 'Position', [p1, 1 - p2 - 0.02, 0.0, -0.02], ...
    'Color', [0, 1, 0]);

if (nCrs > 1)
    psn = zeros(nCrs, 3);
    for t = 1:nCrs
        userdata = get(dataCurs(t).t, 'UserData');
        psn(t, 1:2) = userdata{1};%clicked point (data space)
        psn(t, 3) = userdata{2};%clicked line (handle)
    end
    hstH = unique(psn(:, 3));
    for k = 1:length(hstH)
        jj = (psn(:, 3) == hstH(k));
        psnL = psn(jj, 1:2);
        hL = vertcat(dataCurs(jj).t);
        [~, ii] = sort(psnL(:, 1));
        for z = length(ii):-1:2
            dx = round((psnL(ii(z), 1) - psnL(ii(z - 1), 1)) * 10) / 10;
            dy = round((psnL(ii(z), 2) - psnL(ii(z - 1), 2)) * 10) / 10;
            [~, m] = min(abs(segms(:, 1) - (segms(cs, 1) + xy(1) * timeScalCoeff)));%number of the nearest timestamp
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

PlotAll(handles);%plot with new parameters


% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function Eview_WindowButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to Eview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

evi = handles.Eview;%main form handle
lfpAx = handles.LFP_ax;%LFP axes

x_lim = get(lfpAx, 'XLim');
y_lim = get(lfpAx, 'YLim');
bgnXY = get(lfpAx, 'CurrentPoint');
clkIn = zeros(1, 3);
clkIn(1) = ((bgnXY(1, 1) >= x_lim(1)) && (bgnXY(1, 1) <= x_lim(2)) && ...
            (bgnXY(1, 2) >= y_lim(1)) && (bgnXY(1, 2) <= y_lim(2))); %click point is in axes limits

if ~any(clkIn) %no one axes clicked
    setappdata(evi, 'bgnXY', []);%start point of frame
    if ((bgnXY(1, 1) < x_lim(1)) && ... click was near left edge of LFP axis
        (bgnXY(1, 2) >= y_lim(1)) && (bgnXY(1, 2) <= y_lim(2))) %call for channel group settings
        %number of visible selected channel = round(bgnXY(1, 2))
        TunChannels(round(bgnXY(1, 2)), handles)%call channels tuning dialog
    end
    
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

evi = handles.Eview;%handle of the application window
bgnXY = getappdata(evi, 'bgnXY');%start point of frame
endXY = getappdata(evi, 'endXY');%end point of frame
if (~isempty(bgnXY) && isempty(endXY)) %both point of recangle exist
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

evi = handles.Eview;%handle of the application window
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
    sn = str2double(get(handles.CurrSegm, 'String'));%segments of signal (sweeps)
    
    %rCh = get(handles.FrstCh, 'Value'):get(handles.LastCh, 'Value');%channels to be read and draw
    chTabData = get(handles.ChannSelector, 'Data');
    rCh = 1:hd.nADCNumChannels;
    rCh = rCh(horzcat(chTabData{:, 3}));%selected channels
    
    lfpChlds = get(handles.LFP_ax, 'Children');
    z = 3;
    dCoeff = [];
    while isempty(dCoeff)
        dCoeff = get(lfpChlds(z), 'UserData');
        z = z + 1;
    end
    dCoeff = dCoeff{3};%scale factors
    
    bT.t = segms(sn, 1) + x_lim * dCoeff;%begin-end of currently defined burst (or highlighted area)
    jj = ceil(y_lim(1)):floor(y_lim(2));%highlighted vertical range
    jj = jj((jj > 0) & (jj <= length(rCh)));%selected channels
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

            prm.Fs = zavp.dwnSmplFrq;%discretization frequency for signals to be treated (Hz)
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
        brst = getappdata(evi, 'brst');
        if ~isempty(brst) %bursts exist
            z = length(brst) + 1;%number of the new burst
        else %no one burst yet
            clear brst
            brst = struct('t', [], 'ch', []);%creat bursts data
            z = 1;%number of the new burst
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
            if ((diff(round(bT.t)) > 5) && ~isempty(bT.ch)) %long segment
                %z - number of the new burst
                brst(z).t = round(bT.t);%absolute time (ms = samples)
                brst(z).ch = bT.ch;%channels
                setappdata(evi, 'selectedBrst', z);%selected burst
                setappdata(evi, 'brst', brst);
                set(handles.ShowBrst, 'Enable', 'on')
                set(handles.ShowBrst, 'State', 'on')
                PlotAll(handles);%plot with new parameters
            end
            set(handles.SaveBursts, 'UserData', 1)%set flag of bursts changed
        end
    end %end of (if isequal(get(handles.PwrSpectrWin, 'State'), 'on') %spectrum in window selected)
end


        
function TunChannels(slVisCh, handles)
%call dialog for channels tuning
%
%slVisCh - number of selected visible channel
%handles - structure with handles and user data (see GUIDATA)
%
evi = handles.Eview;%handle of the application window
chanSettingGrp = getappdata(evi, 'chanSettingGrp');%load channels group setting
chTabData = get(handles.ChannSelector, 'Data');%read current table content

numSelCh = find(horzcat(chTabData{:, 3}));%number of selected channels
slctCh = numSelCh(slVisCh);%original number of selected channel

%= creat dialog window =%
tunFig = figure('Name', 'Tune channels view', 'Units', 'normalized', 'MenuBar','none', ...
    'NumberTitle', 'off', 'Position',[0.40546875 0.6455078125 0.3421875 0.1357421875], ...
    'Resize', 'off', 'Tag', 'TunChanFig', 'Visible', 'off', 'WindowStyle', 'modal', 'Color', [0.831, 0.816, 0.784]);

tfh(1:21) = struct('h', []);%structure with tune figure control handles
tfh(1).h = uicontrol('Parent', tunFig, 'Units', 'normalized', 'Tag', 'ApplySets', 'Style', 'pushbutton', ...
    'Position', [0.86986301369863 0.0359712230215827 0.107305936073059 0.158273381294964], ...
    'String', 'OK', 'Callback', {@ApplyNewChanSets, handles}, 'BackgroundColor', [0.831, 0.816, 0.784]);

tfh(2).h = uicontrol('Parent', tunFig, 'Units', 'normalized', 'String', num2str(chanSettingGrp(slctCh).lowCutFrq), ...
    'Style', 'edit', 'Tag', 'LowCutFrq', 'BackgroundColor', [1, 1, 1], ...
    'Position', [0.0228310502283105 0.525179856115108 0.136986301369863 0.143884892086331]);

tfh(3).h = uicontrol('Parent', tunFig, 'Units','normalized', 'String', num2str(chanSettingGrp(slctCh).highCutFrq), ...
    'Style', 'edit', 'Tag', 'HighCutFrq', 'BackgroundColor', [1, 1, 1], ...
    'Position', [0.271689497716894 0.525179856115108 0.136986301369863 0.143884892086331]);

tfh(4).h = uicontrol('Parent', tunFig, 'Units', 'normalized', 'String', 'Hz', 'Style', 'text', ...
    'Tag', 'text2', 'BackgroundColor', [0.831, 0.816, 0.784], ...
    'Position', [0.17351598173516 0.546762589928059 0.045662100456621 0.100719424460432]);

tfh(5).h = uicontrol('Parent', tunFig, 'Units','normalized', 'String', 'Hz', 'Style', 'text', ...
    'Tag', 'text3','BackgroundColor', [0.831, 0.816, 0.784], ... 
    'Position', [0.422374429223744 0.546762589928059 0.045662100456621 0.100719424460432]);

tfh(6).h = uicontrol('Parent', tunFig, 'Units','normalized', 'String', 'Lowcut filter', 'Style', 'checkbox', ...
    'Position', [0.0228310502283105 0.719424460431657 0.212328767123288 0.115107913669065], 'FontSize', 9, ...
    'Tag', 'LowCutFilter', 'Value', chanSettingGrp(slctCh).lowCut, 'BackgroundColor', [0.831, 0.816, 0.784]);

tfh(7).h = uicontrol('Parent', tunFig, 'Units', 'normalized', 'String', 'Highcut filter', 'Style', 'checkbox', ...
    'Position', [0.271689497716894 0.719424460431657 0.212328767123288 0.115107913669065], 'FontSize', 9, ...
    'Tag', 'HighCutFilter', 'Value', chanSettingGrp(slctCh).highCut, 'BackgroundColor', [0.831, 0.816, 0.784]);

tfh(8).h = uicontrol('Parent', tunFig, 'Units', 'normalized', 'FontSize', 10, 'Style', 'text', 'Tag','ChGrpName', ...
    'Position', [0.0114155251141553 0.856115107913669 0.977168949771689 0.122302158273381], 'UserData', slctCh, ...
    'HorizontalAlignment', 'left', 'String', ['channel group: ', num2str(chanSettingGrp(slctCh).grpName)], 'BackgroundColor', [0.9, 0.9, 0.9]);

tfh(9).h = uicontrol('Parent', tunFig, 'Units', 'normalized', 'String', num2str(chanSettingGrp(slctCh).scale), ...
    'Position', [0.289497716894977 0.16546762589928066 0.157534246575342 0.143884892086331], ...
    'Style', 'edit', 'Tag', 'SgnlScale', 'BackgroundColor', [1, 1, 1]);

tfh(10).h = uicontrol('Parent', tunFig, 'Units', 'normalized', 'Style', 'text', 'FontSize', 9, ...
    'Position', [0.449315068493151 0.18705035971223044 0.045662100456621 0.100719424460432], ...
    'String', 'mV', 'Tag', 'text5', 'BackgroundColor', [0.831, 0.816, 0.784]);

tfh(11).h = uicontrol('Parent', tunFig, 'Units', 'normalized', 'Style', 'text', 'FontSize', 10, ...
    'Position', [0.264383561643836 0.3093525179856118 0.205479452054795 0.122302158273381], ...
    'String', 'vertical scale', 'Tag', 'text6', 'BackgroundColor', [0.831, 0.816, 0.784]);

tfh(12).h = uicontrol('Parent', tunFig, 'Units', 'normalized', 'String', 'Comm.ref.', 'Style', 'checkbox', ...
    'Position', [0.801369863013698 0.446043165467627 0.182648401826484 0.14], 'FontSize', 9, ...
    'Tag', 'CommRef', 'Value', chanSettingGrp(slctCh).comRef, 'BackgroundColor', [0.831, 0.816, 0.784]);

tfh(13).h = uicontrol('Parent', tunFig, 'Units', 'normalized', 'String', 'CSD map', 'Style', 'checkbox',...
    'Position', [0.557077625570776 0.18705035971223 0.178082191780822 0.136690647482014], ...
    'Callback', {@CSDmaporMUAmap, handles}, 'FontSize', 9, ...
    'Tag', 'PlotCSDmap', 'Value', chanSettingGrp(slctCh).pltCSDmap, 'BackgroundColor', [0.831, 0.816, 0.784]);

tfh(14).h = uicontrol('Parent', tunFig, 'Units', 'normalized', 'String', 'DC shift', 'Style', 'checkbox', ...
    'Position', [0.801369863013699 0.697841726618705 0.166666666666667 0.136690647482014], 'FontSize', 9, ...
    'Tag', 'DCshift', 'Value', chanSettingGrp(slctCh).dcShift, 'BackgroundColor', [0.831, 0.816, 0.784]);

tfh(15).h = uicontrol('Parent', tunFig, 'Units', 'normalized', 'String', 'raw signal', 'Style', 'checkbox',...
    'Position',[0.0319634703196347 0.129496402877698 0.198630136986301 0.165467625899281], 'FontSize', 9, ...
    'Tag', 'RawSignal', 'Value', chanSettingGrp(slctCh).rawSignal, 'BackgroundColor', [0.831, 0.816, 0.784]);

tfh(16).h = uicontrol('Parent', tunFig, 'Units', 'normalized', 'String', 'Spikes map', 'Style', 'checkbox',...
    'Position', [0.557077625570776 0.0359712230215828 0.21689497716895 0.136690647482014], ...
    'Callback', {@CSDmaporMUAmap, handles}, 'FontSize', 9, ...
    'Tag', 'PlotSpkMap', 'Value', chanSettingGrp(slctCh).pltSpkMap, 'BackgroundColor', [0.831, 0.816, 0.784]);

tfh(17).h = uicontrol('Parent', tunFig, 'Units', 'normalized', 'Style', 'pushbutton', ...
    'Position', [0.287671232876712 0.0287769784172662 0.0662100456621005 0.129496402877698],...
    'Callback', {@DownsizeAmplScale, handles}, ...
    'String', '<', 'Tag', 'DownsizeAmpl', 'BackgroundColor', [0.831, 0.816, 0.784]);

tfh(18).h = uicontrol('Parent', tunFig, 'Units', 'normalized', 'Style', 'pushbutton', ...
    'Position', [0.378995433789954 0.0287769784172662 0.0662100456621005 0.13], ...
    'Callback', {@UpsizeAmplScale, handles}, ...
    'String', '>', 'Tag', 'UpsizeAmpl', 'BackgroundColor', [0.831, 0.816, 0.784]);

tfh(19).h = uicontrol('Parent', tunFig, 'Units', 'normalized', 'String', 'MUA', 'Style', 'checkbox', ...
    'Position', [0.557077625570776 0.697841726618705 0.125570776255708 0.136690647482014], 'FontSize', 9, ...
    ...'Callback', {@MUAorMUAfreq, handles},...
    'Tag', 'PlotMUA', 'Value', chanSettingGrp(slctCh).pltMUA, 'BackgroundColor', [0.831, 0.816, 0.784]);

tfh(20).h = uicontrol('Parent', tunFig, 'Units', 'normalized', 'String', 'LFP', 'Style','checkbox', ...
    'Position', [0.557077625570776 0.410071942446044 0.125570776255708 0.136690647482014], 'FontSize', 9, ...
    'Tag', 'PlotLFP', 'Value', chanSettingGrp(slctCh).pltLFP, 'BackgroundColor', [0.831, 0.816, 0.784]);

tfh(21).h = uicontrol('Parent', tunFig, 'Units', 'normalized', 'String', 'MUA freq', 'Style', 'checkbox', ...
    'Position', [0.557077625570776 0.553956834532374 0.182648401826484 0.136690647482014], 'FontSize', 9, ...
    ...'Callback', {@MUAorMUAfreq, handles}, ...
    'Tag', 'PlotMUAfreq', 'Value', chanSettingGrp(slctCh).pltMUAfreq, 'BackgroundColor', [0.831, 0.816, 0.784]);
%= end of dialog window creation =%

%transfer dialog-control handles
setappdata(evi, 'tunChFigH', tunFig);%save tune-channel window handle
setappdata(evi, 'tunFigHndls', tfh);%save tune figure handles

set(tunFig, 'Visible', 'on')%show dialog



% --- accept setting of channels group and apply --- %
function ApplyNewChanSets(hObject, eventdata, handles)
% hObject    handle to Merge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
%accept setting of channels group and apply
evi = handles.Eview;%handle of the application window
tunFig = getappdata(evi, 'tunChFigH');%load tune-channel window handle
tfh = getappdata(evi, 'tunFigHndls');%load handles of contorl of tune-channel window
chanSettingGrp = getappdata(evi, 'chanSettingGrp');%load channels group setting
slctCh = get(tfh(8).h, 'UserData');%original number of selected channel

for t = 1:length(chanSettingGrp) %run over channels (looking for the group)
    if (chanSettingGrp(t).grpName == chanSettingGrp(slctCh).grpName) %same group
        chanSettingGrp(t).lowCut = get(tfh(6).h, 'Value');%if do lowcut filtering
        chanSettingGrp(t).lowCutFrq = abs(str2double(get(tfh(2).h, 'String')));%lowcut frequency
        chanSettingGrp(t).highCut = get(tfh(7).h, 'Value');%if do highcut filtering
        chanSettingGrp(t).highCutFrq = abs(str2double(get(tfh(3).h, 'String')));%highcut frequency
        chanSettingGrp(t).rawSignal = get(tfh(15).h, 'Value');%if plot raw signal
        chanSettingGrp(t).comRef = get(tfh(12).h, 'Value');%if do common reference subtraction
        chanSettingGrp(t).dcShift = get(tfh(14).h, 'Value');%if compensate DC
        chanSettingGrp(t).pltCSDmap = get(tfh(13).h, 'Value');%if plot CSD map
        chanSettingGrp(t).pltSpkMap = get(tfh(16).h, 'Value');%if plot spikes map
        chanSettingGrp(t).pltLFP = get(tfh(20).h, 'Value');%plot LFP traces
        chanSettingGrp(t).pltMUA = get(tfh(19).h, 'Value');%plot MUA lines
        chanSettingGrp(t).pltMUAfreq = get(tfh(21).h, 'Value');%plot trace of spikes densities
        chanSettingGrp(t).scale = abs(str2double(get(tfh(9).h, 'String')));%verticale scale (mV, amplitude)
    end
end

delete(tunFig)
setappdata(evi, 'chanSettingGrp', chanSettingGrp);%save channels group setting
PlotAll(handles);%plot with new parameters



function CSDmaporMUAmap(hObject, eventdata, handles)
% hObject    handle to Merge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
%choose one of map: CSD or MUA
if get(hObject, 'Value') %one of map was choosed
    evi = handles.Eview;%handle of the application window
    tfh = getappdata(evi, 'tunFigHndls');%load handles of contorl of tune-channel window
    if isequal(get(hObject, 'Tag'), 'PlotCSDmap')
        set(tfh(16).h, 'Value', 0)%tfh(16).h -> PlotSpkMap
    elseif isequal(get(hObject, 'Tag'), 'PlotSpkMap')
        set(tfh(13).h, 'Value', 0)%tfh(13).h -> PlotCSDmap
    end
end



function MUAorMUAfreq(hObject, eventdata, handles)
% hObject    handle to Merge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
%choose one of line type: MUA (verticale lines) or MUA frequency (timecourse)
if get(hObject, 'Value') %one of map was choosed
    evi = handles.Eview;%handle of the application window
    tfh = getappdata(evi, 'tunFigHndls');%load handles of contorl of tune-channel window
    if isequal(get(hObject, 'Tag'), 'PlotMUA')
        set(tfh(21).h, 'Value', 0)%tfh(21).h -> PlotMUAfreq
    elseif isequal(get(hObject, 'Tag'), 'PlotMUAfreq')
        set(tfh(19).h, 'Value', 0)%tfh(19).h -> PlotMUA
    end
end



% --- change LFP scale (downsize) --- %
function DownsizeAmplScale(hObject, eventdata, handles)
% hObject    handle to Merge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
evi = handles.Eview;%handle of the application window
tfh = getappdata(evi, 'tunFigHndls');%load handles of contorl of tune-channel window
vScl = abs(str2double(get(tfh(9).h, 'String'))) / 2;%new scale (downsized)
if (vScl >= 0.1)
    vScl = round(10 * vScl) / 10;
elseif ((vScl < 0.1) && (vScl >= 0.01))
    vScl = round(100 * vScl) / 100;
else
    vScl = round(1000 * vScl) / 1000;
end
if (vScl <= 0)
    vScl = 0.001;
end
set(tfh(9).h, 'String', num2str(vScl))


% --- change LFP scale (upsize) --- %
function UpsizeAmplScale(hObject, eventdata, handles)
% hObject    handle to Merge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
evi = handles.Eview;%handle of the application window
tfh = getappdata(evi, 'tunFigHndls');%load handles of contorl of tune-channel window
vScl = abs(str2double(get(tfh(9).h, 'String'))) * 2;%new scale (upsized)
if (vScl >= 0.1)
    vScl = round(10 * vScl) / 10;
elseif ((vScl < 0.1) && (vScl >= 0.01))
    vScl = round(100 * vScl) / 100;
else
    vScl = round(1000 * vScl) / 1000;
end
if (vScl <= 0)
    vScl = 0.001;
end
set(tfh(9).h, 'String', num2str(vScl))


% --------------------------------------------------------------------
function DeleteBrst_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to DeleteBrst (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
evi = handles.Eview;%handle of the application window
brst = getappdata(evi, 'brst');
if ~isempty(brst) %bursts exist
    slB = getappdata(evi, 'selectedBrst');%selected burst
    if ~isempty(slB) %any burst selected
        brst(slB) = [];%delete selected burst
        if isempty(brst)
            set(handles.ShowBrst, 'Enable', 'off')
        end
        setappdata(evi, 'brst', brst);
        setappdata(evi, 'selectedBrst', []);
        set(handles.SaveBursts, 'UserData', 1)%set flag of bursts changed
        PlotAll(handles);%plot with new parameters
    end
end


% --- Executes when user attempts to close Eview.
function Eview_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to Eview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

evi = handles.Eview;%main form handle
if (~isempty(getappdata(evi, 'brst')) && get(handles.SaveBursts, 'UserData')) %flag of bursts changed
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
brst = getappdata(evi, 'brst');
if isempty(brst) %bursts exist
    pth = getappdata(evi, 'pth');%directory
    flNm = getappdata(evi, 'flNm');%filename
    
    lfp = getappdata(evi, 'lfp');
    %lfp = double(lfp);
    spks = getappdata(evi, 'spks');
    %lfpVar = getappdata(evi, 'lfpVar');
    hd = getappdata(evi, 'hd');
    zavp = getappdata(evi, 'zavp');
    prm(1:length(brst)) = struct('p', []);
    jj = vertcat(brst(:).t);
    [~, ii] = sort(jj(:, 1));
    brst = brst(ii);%sorted bursts
    
    filtFlag = false(hd.nADCNumChannels, 1);%true if lfp was fitered
    for nch = 1:hd.nADCNumChannels
        try
            cscHd = Nlx2MatCSC([pth, flNm, '\', hd.recChNames{nch, 1}, '.ncs'], [0 0 0 0 0], 1, 1, []);%header (lfp)
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
    fprintf(fid, 'Window length (ms)\t%d\n', str2double(get(handles.AftStim, 'String')) + str2double(get(handles.BefStim, 'String')));%time window length
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
    fprintf(fid, '#\tms\tms\t#\tms\tms\tuV\tms\tms\tspikes\n');
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
set(handles.SaveBursts, 'UserData', 0)%set flag of bursts saved
close(hWt);%destroy waitbar



function HistPrecision_Callback(hObject, eventdata, handles)
% hObject    handle to HistPrecision (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

PlotAll(handles);%plot with new parameters


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
[imgFlNm, pth] = uiputfile({'*.emf', 'EMF'; '*.jpg', 'JPG'}, 'Save as', [pth, flNm]);%save image dialog
z = find(imgFlNm == '.', 1, 'last');

if isequal(imgFlNm((z + 1):end), 'emf') %EMF
    saveas(evi, [pth, imgFlNm], imgFlNm((z + 1):end))%save image
elseif isequal(imgFlNm((z + 1):end), 'jpg') %JPG
    imgM = getframe(evi);%get objects of picture
    imwrite(imgM.cdata, [pth, imgFlNm], imgFlNm((z + 1):end))%save image
end



% --- Executes on button press in SelectChann.
function SelectChann_Callback(hObject, eventdata, handles)
% hObject    handle to SelectChann (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
%channels visibility settings
chanSelecH = handles.ChannSelector;
if isequal(get(chanSelecH, 'Visible'), 'on')
    set(chanSelecH, 'Visible', 'off')
    set(handles.SelectAllCh, 'Visible', 'off')
    set(handles.InvertSelection, 'Visible', 'off')
    PlotAll(handles);%plot with new parameters
else
    set(chanSelecH, 'Visible', 'on')
    set(handles.SelectAllCh, 'Visible', 'on')
    set(handles.InvertSelection, 'Visible', 'on')
end



function MUAScale_Callback(hObject, eventdata, handles)
% hObject    handle to MUAScale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

PlotAll(handles);%plot with new parameters


% --- Executes on button press in SelectAllCh.
function SelectAllCh_Callback(hObject, eventdata, handles)
% hObject    handle to SelectAllCh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

chanSelecH = handles.ChannSelector;%table of channels selections
chTabData = get(chanSelecH, 'Data');%data of the table
for t = 1:size(chTabData, 1) %run over rows of table (channels)
    chTabData{t, 3} = true;%select all channels
end
set(chanSelecH, 'Data', chTabData)%set new table data



% --- Executes on button press in InvertSelection.
function InvertSelection_Callback(hObject, eventdata, handles)
% hObject    handle to InvertSelection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

chanSelecH = handles.ChannSelector;%table of channels selections
chTabData = get(chanSelecH, 'Data');%data of the table
for t = 1:size(chTabData, 1) %run over rows of table (channels)
    chTabData{t, 3} = ~chTabData{t, 2};%invert the selection
end
set(chanSelecH, 'Data', chTabData)%set new table data


% --- Executes on button press in WinAvr.
function WinAvr_Callback(hObject, eventdata, handles)
% hObject    handle to WinAvr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(hObject, 'Value')
    set(handles.Average, 'Value', 1);
    PlotAll(handles);%plot with new parameters
end



% --------------------------------------------------------------------
function SaveScreenMat_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to SaveScreenMat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%export data of displayed traces to mat-file
lfpTrcs = findobj(get(handles.LFP_ax, 'Children'), 'Type', 'line');
lfpTrcs = findobj(lfpTrcs, 'Color', [0, 0, 0]);

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
    k = k + 1;
end
tm_p = tm_p * dCoeffP{3};%lfp time vector in milliseconds
if spkshown %spikes visible
    tm_s = tm_s * dCoeffS{3};%spikes time vector in milliseconds
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




% --- Executes when selected cell(s) is changed in ChannSelector.
function ChannSelector_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to ChannSelector (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)

ii = unique(eventdata.Indices(:, 1));%indices of selected rows
set(handles.ChannSelector, 'UserData', eventdata.Indices)%save current selection
if ((length(ii) > 1) && all(diff(ii) == 1) && all(eventdata.Indices(:, 2) == 3)) %change channels selection
    chTabData = get(handles.ChannSelector, 'Data');%read current table content
    jj = horzcat(chTabData{ii, 3});%selection falgs
    trg = ~(sum(jj) >= sum(~jj));%target value (prevailing selection)
    for z = ii' %run over selected raws
        chTabData{z, 3} = trg;%change value of selected cells (set prevailing value)
    end
    set(handles.ChannSelector, 'Data', chTabData)%set table data
    
%     jscrollpane = javaObjectEDT(get(handles.ChannSelector, 'UserData'));%java-handle of scrollpane
%     viewport = javaObjectEDT(jscrollpane.getViewport);%jave-handle
%     jtable = javaObjectEDT(viewport.getView);%jave-handle of table content
%     jtable.scrollRowToVisible(ii(end));%make selected rows visible again
end

if ((length(ii) > 1) && all(diff(ii) == 1) && all(eventdata.Indices(:, 2) == 4)) %change channels groups
    evi = handles.Eview;%handle of the application window
    chTabData = get(handles.ChannSelector, 'Data');%read current table content
    chanSettingGrp = getappdata(evi, 'chanSettingGrp');%load channels group setting
    
    ngc = vertcat(chTabData{ii, 4});%names of groups of selected channels
    [unqngc, ~, jj] = unique(ngc);%positions of unique names
    if (length(unique(jj)) == 1) %separat selected channels into new group
        m = min(setdiff(1:size(chTabData, 1), horzcat(chTabData{:, 4})));%number of new unique group to be added
        t = ii(1);%number of the first selected channel (template settings)
    else %more than one group are merged
        m = unqngc(mode(jj));%name of group for all selected channels (group prevailing by number of channels)
        for z = 1:length(chanSettingGrp) %run over channel
            if (chanSettingGrp(z).grpName == m)
                t = z;%number of the first channel in group (template settings)
                break;%out of (for z)
            end
        end
    end
    for z = ii' %run over selected rows (channels)
        chanSettingGrp(z).grpName = m;%name of group for all selected channels
        chTabData{z, 4} = chanSettingGrp(z).grpName;%copy group name
        chanSettingGrp(z) = chanSettingGrp(t);%copy setting from template
    end
    set(handles.ChannSelector, 'Data', chTabData)%set table data
    setappdata(evi, 'chanSettingGrp', chanSettingGrp);%save channels group setting
end


function SeparateChannel(hObject, eventdata, handles)
%creat individual group for the selected channel
%
evi = handles.Eview;%handle of the application window
chTabData = get(handles.ChannSelector, 'Data');%read current table content
chanSettingGrp = getappdata(evi, 'chanSettingGrp');%load channels group setting
ii = get(handles.ChannSelector, 'UserData');%load last selection
ii = unique(ii(:, 1));%indices of selected rows
m = min(setdiff(1:size(chTabData, 1), horzcat(chTabData{:, 4})));%number of new unique group to be add

for z = ii' %run over selected rows (channels)
    chanSettingGrp(z).grpName = m;%name of group for all selected channels
    chTabData{z, 4} = chanSettingGrp(z).grpName;%copy group name
end
set(handles.ChannSelector, 'Data', chTabData)%set table data
setappdata(evi, 'chanSettingGrp', chanSettingGrp);%save channels group setting
    


function CSDcMapMax_Callback(hObject, eventdata, handles)
% hObject    handle to CSDcMapMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%change color axis
set(handles.CSDcMapMin, 'String', ['-', get(handles.CSDcMapMax, 'String')]) %set minimum automatically
caxis(handles.LFP_ax, [str2double(get(handles.CSDcMapMin, 'String')) / 1e3, str2double(get(handles.CSDcMapMax, 'String')) / 1e3])%colorbar range



function SpkCMap_Callback(hObject, eventdata, handles)
% hObject    handle to SpkCMapMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%change color axis
caxis(handles.LFP_ax, [str2double(get(handles.SpkCMapMin, 'String')), str2double(get(handles.SpkCMapMax, 'String'))])%colorbar range
