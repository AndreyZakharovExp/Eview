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

% Last Modified by GUIDE v2.5 21-May-2021 18:08:53

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

% set(evi, 'WindowButtonDownFcn', {@Eview_WindowButtonDownFcn, handles})%reaction on button down
% set(evi, 'WindowButtonMotionFcn', {@Eview_WindowButtonMotionFcn, handles})%reaction on button motion
% set(evi, 'WindowButtonUpFcn', {@Eview_WindowButtonUpFcn, handles})%reaction on button up

h = uicontextmenu;%contex
uimenu(h, 'Label', 'Set unique', 'Callback', {@SeparateChannel, handles});
set(handles.ChannSelector, 'UIContextMenu', h);%associate context menu with channel-seletor table

manDetPrmH = handles.ManDetPrm;%table of manual detection parameters
set(manDetPrmH, 'ColumnName', {'PrmName', 'Value', 'Description'})
set(manDetPrmH, 'ColumnEditable', [false, true, false])
set(manDetPrmH, 'ColumnFormat', {'char', 'numeric', 'char'})
set(manDetPrmH, 'ColumnWidth', {80, 40, 220})
set(manDetPrmH, 'Data', {'Channel', 1, 'main channel for manual detection'; ... main channel for manual detection
                         'Polarity', 1, '1-min, 0-max'; ... extremum type
                         'Search range', 15, 'range for local extremum search'; ... extremum search range
                         'SubChannel', 2, 'subchannel for additional analysis'; ... minor channel 
                         'SubPolarity', 1, 'polarity on subchannel, 1-min, 0-max' ... polarity minor channel 
    })%set manual detection parameters

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
    return;%exit from OpenData_ClickedCallback
end

z = find(flNm == '.', 1, 'last');
if ~isequal(flNm(z:end), '.mat')
    %set(nts, 'String', 'unknown file type')
    return;%exit from OpenData_ClickedCallback
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
%= end of (load events from xls) =%

%= detected bursts (any segments of recordation) =%
brst = [];%clear previous bursts
set(handles.ShowBrst, 'Enable', 'off')
if exist([pth, flNm, '_brst.mat'], 'file')
    load([pth, flNm, '_brst.mat'])
    if ~isstruct(brst) %wrong format
        jj = brst;%copy
        clear brst %renew variable
        brst(1:size(jj, 1)) = struct('t', [], 'ch', []);
        for t = 1:size(jj, 1) %run over bursts
            brst(t).t = jj(t, 1:2);%begin-end time moments of bursts
            brst(t).ch = 1:hd.nADCNumChannels;%channels of bursts
        end
    end
    set(handles.ShowBrst, 'Enable', 'on')%burst can be shown
end
set(handles.SaveBursts, 'UserData', 0)%set flag of bursts saved
%= end of (detected bursts (any segments of recordation)) =%

%= fill table of channel selector =%
chTabData = get(handles.ChannSelector, 'Data');%channel selection table
if isempty(horzcat(chTabData{:, 1})) %first initiation
    chTabData = cell(hd.nADCNumChannels, 2);%create channels-table data
    
    %creat channel groups with default settings 
    chanSettingGrp(1:hd.nADCNumChannels) = struct('grpName', [], ... name of channels group
                                                  'lowCut', false, ... if od lowcut filtering
                                                  'lowCutFrq', 0.1, ... lowcut frequency
                                                  'highCut', false, ... if do highcut filtering
                                                  'highCutFrq', 100, ... highcut frequency
                                                  'recovDC', false, ... if recovery DC signal from psedoDC
                                                  'rawSignal', false, ... if raw signal
                                                  'comRef', false, ... if common reference
                                                  'dcShift', false, ... if compensate DC
                                                  'invertLFP', false, ... if invert LFP
                                                  'pltCSDmap', false, ... if plot CSD map
                                                  'pltSpkMap', false, ... if plot spikes frequency map
                                                  'pltLFP', true, ... if plot LFP traces
                                                  'pltMUA', true, ... if plot MUA bars
                                                  'pltMUAfreq', false, ... if plot MUA-frequency traces
                                                  'lfpScale', 1, ... LFP-scale (mV)
                                                  'muaTrsld', 5, ... MUA threshold (std)
                                                  'muaBin', 5, ... bin step for MUA frequency (ms)
                                                  'precisF', 1, ... precision factor (for MUA frequency calculation)
                                                  'muaScale', 1, ... MUA-scale (spikes/ms)
                                                  'realsernum', 0 ... actual serial number of the channel on a graph (number of line)
                                                  );
    for ch = 1:hd.nADCNumChannels %run over channels
        chTabData{ch, 3} = false;%channel deselected
        chanSettingGrp(ch).grpName = ch;%initially each channel is in its own group
        chTabData{ch, 4} = chanSettingGrp(ch).grpName;%copy group name
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
                                                  'lowCut', false, ... if od lowcut filtering
                                                  'lowCutFrq', 0.1, ... lowcut frequency
                                                  'highCut', false, ... if do highcut filtering
                                                  'highCutFrq', 100, ... highcut frequency
                                                  'recovDC', false, ... if recovery DC signal from psedoDC
                                                  'rawSignal', false, ... if raw signal
                                                  'comRef', false, ... if common reference
                                                  'dcShift', false, ... if compensate DC
                                                  'invertLFP', false, ... if invert LFP
                                                  'pltCSDmap', false, ... if plot CSD map
                                                  'pltSpkMap', false, ... if plot spikes frequency map
                                                  'pltLFP', true, ... if plot LFP traces
                                                  'pltMUA', true, ... if plot MUA bars
                                                  'pltMUAfreq', false, ... if plot MUA-frequency traces
                                                  'lfpScale', 1, ... LFP-scale (mV)
                                                  'muaTrsld', 5, ... MUA threshold (std)
                                                  'muaBin', 5, ... bin step for MUA frequency (ms)
                                                  'precisF', 1, ... precision factor (for MUA frequency calculation)
                                                  'muaScale', 1, ... MUA-scale (spikes/ms)
                                                  'realsernum', 0 ... actual serial number of the channel on a graph (number of line)
                                                  );%add structures
        m = max(horzcat(chanSettingGrp(1:z).grpName)) + 1;%new groups for each new channel
        for ch = (z + 1):hd.nADCNumChannels %run over channels
            chTabData{ch, 3} = false;%channel deselected
            chanSettingGrp(ch).grpName = m;%each new channel is in its own group
            chTabData{ch, 4} = chanSettingGrp(ch).grpName;%copy group name
            m = m + 1;%new group for each new channel
        end
    end
    for ch = 1:hd.nADCNumChannels %run over all channels
        chanSettingGrp(ch).rawSignal = false;%no raw signal draw
    end
end

%rewrite names of channels
if isfield(hd, 'recChNames') %channel names exist
    for ch = 1:hd.nADCNumChannels %run over channels
        chTabData{ch, 1} = hd.recChNames{ch, 1};%name of channel
        chTabData{ch, 2} = '';%altername of channel
    end
    if (size(hd.recChNames, 2) > 1) %altername was signed
        for ch = 1:hd.nADCNumChannels %run over channels
            chTabData{ch, 2} = hd.recChNames{ch, 2};%name of channel
        end
    end
else %channel names not exist
    hd.recChNames = cell(hd.nADCNumChannels, 1);
    for ch = 1:hd.nADCNumChannels %run over channels
        chTabData{ch, 1} = num2str(ch);%name of channel
    end
end
if ~any(vertcat(chTabData{:, 3})) %no one channel selection
    chTabData{1, 3} = true;%select only first channel
end
set(handles.ChannSelector, 'Data', chTabData)%set table data
%= end of (fill table of channel selector) =%

%= LFP/spikes =%
lfp = double(lfp);
spksOrig = spks;%copy of spikes
for sw = 1:hd.lActualEpisodes
    for ch = 1:hd.nADCNumChannels
        prg = chanSettingGrp(ch).muaTrsld;%threshold for the channel
        if isfield(spksOrig, 's')
            spksOrig(ch, sw).tStamp = spksOrig(ch, sw).s(:);%
            spks(ch, sw).tStamp = spks(ch, sw).s(:);%
        end
        if isfield(spksOrig, 'ampl')
            ii = double(spksOrig(ch, sw).ampl) <= (-lfpVar(ch) * prg);
            spks(ch, sw).tStamp = spksOrig(ch, sw).tStamp(ii);
            spks(ch, sw).ampl = spksOrig(ch, sw).ampl(ii);
        end
        if isempty(spks(ch, sw).tStamp)
            spks(ch, sw).tStamp = zeros(0, 1);
            spks(ch, sw).ampl = zeros(0, 1);
        end
    end
end
%= end of (LFP/spikes) =%

if iscell(zavp.file)
    lfp = RemovNaN(lfp);%exclude NaN and Inf with interpolation
end
lfpOrig = lfp;%original LFP

%= LFP filtration =%
for ch = 1:length(chanSettingGrp) %run over channels (looking for the group)
    if (chanSettingGrp(ch).recovDC) %recovery DC from pseudoDC
        lfp(:, ch, :) = RecovDCfromPseudo(lfp(:, ch, :), zavp.dwnSmplFrq);%recovery
    end
    if (chanSettingGrp(ch).lowCut && (chanSettingGrp(ch).lowCutFrq < (zavp.dwnSmplFrq / 2))) %lowcut filtration
        if (chanSettingGrp(ch).lowCutFrq < 2) %small cutting frequency
            lfp(:, ch, :) = ZavRCfilt(lfp(:, ch, :), chanSettingGrp(ch).lowCutFrq, zavp.dwnSmplFrq, 'high');%RC-filter, lowcut
        else %quite large cutting frequency
            lfp(:, ch, :) = ZavFilter(lfp(:, ch, :), zavp.dwnSmplFrq, 'high', chanSettingGrp(ch).lowCutFrq, 2);%lowcut
        end
    end
    if (chanSettingGrp(ch).highCut && (chanSettingGrp(ch).highCutFrq < (zavp.dwnSmplFrq / 2))) %highcut filtration
        if (chanSettingGrp(ch).highCutFrq < 2) %small cutting frequency
            lfp(:, ch, :) = ZavRCfilt(lfp(:, ch, :), chanSettingGrp(ch).highCutFrq, zavp.dwnSmplFrq, 'low');%RC-filter, lowcut
        else
            lfp(:, ch, :) = ZavFilter(lfp(:, ch, :), zavp.dwnSmplFrq, 'low', chanSettingGrp(ch).highCutFrq, 2);%highcut
        end
    end
    
    if (chanSettingGrp(ch).invertLFP) %if invert LFP
        lfp(:, ch, :) = -1 * lfp(:, ch, :);%invert LFP
    end
end
%= end of (LFP filtration) =%

%= line handles =%
lfp_ax = handles.LFP_ax;%LFP axis handle
cla(lfp_ax)%clear LFP-axis
lfpLn = zeros(hd.nADCNumChannels, 1);%LFP line handles
spkBr = zeros(hd.nADCNumChannels, 1);%spike bar handles
spkDnLn = zeros(hd.nADCNumChannels, 1);%spike density line handles
clrMapH = zeros(hd.nADCNumChannels, 1);%colormap handles (CSD or MUA)
for ch = 1:hd.nADCNumChannels %run over channels
    subplot(lfp_ax);%set lfp_ax as a parent for next plot
    clrMapH(ch) = imagesc(0, 0, 1);%colormap handles (CSD or MUA)
    set(clrMapH(ch), 'Visible', 'off', 'Tag', 'clrMap')
end
for ch = 1:hd.nADCNumChannels %run over channels
    lfpLn(ch) = plot(lfp_ax, NaN, NaN, 'k', 'LineWidth', 1);%LFP line handles
    set(lfpLn(ch), 'ButtonDownFcn', {@Line_ButtonDownFcn, handles}, 'Visible', 'off', 'Tag', 'lfpLn')%reaction on mouse-click
    spkBr(ch) = plot(lfp_ax, NaN, NaN, 'r', 'LineWidth', 1.5);%spike bar handles
    set(spkBr(ch), 'Visible', 'off')
    spkDnLn(ch) = plot(lfp_ax, NaN, NaN, 'LineWidth', 1);%spike density line handles
    set(spkDnLn(ch), 'ButtonDownFcn', {@Line_ButtonDownFcn, handles}, 'Visible', 'off')%reaction on mouse-click
end
%= end of (line handles) =%

%= save application data =%
setappdata(evi, 'pth', pth) %directory
setappdata(evi, 'flNm', flNm) %filename
setappdata(evi, 'lfp', lfp) %LFP
setappdata(evi, 'lfpOrig', lfpOrig) %original LFP
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
setappdata(evi, 'manlDet', []);%manually detected events
setappdata(evi, 'lfpLn', lfpLn);%LFP line handles
setappdata(evi, 'spkLn', spkBr);%spke bar handles
setappdata(evi, 'spkDnLn', spkDnLn);%spike density line handles
setappdata(evi, 'clrMapH', clrMapH);%colormap handles (CSD or MUA)
setappdata(evi, 'dltGrphObj', []);%deletable objects

set(handles.CutManual, 'Value', 0);%switch off manually detected events
set(handles.CutManual, 'String', 'mdet ');
%= end of (save application data) =%

extSegms = zeros(numel(vertcat(zavp.realStim(:).r)), 4);%preallocation of memory for stimuli moments
z = 1;%through counter of stimuli
for sw = 1:hd.lActualEpisodes %run over sweeps
    for t = 1:numel(zavp.realStim(sw).r) %run over synchro-events
        extSegms(z, 1) = zavp.realStim(sw).r(t) / zavp.rarStep(1);%position of synchro-event (ms from record begin (from 0))
        %extSegms(z, 1) = segms(z, 1) + sepOnsetPeak(15, sw).r(z, 1); disp(' shifted start point ') %shift to SEP onset
        extSegms(z, 2) = zavp.stimCh;%number of channel where syncro-event was detected
        extSegms(z, 3) = sw;%number of sweep where syncro-event was detected
        extSegms(z, 4) = t;%number of stimulus (in matrix zavp.realStim) or number of trough in brst matrix
        z = z + 1;%through counter of stimuli
    end
end
extSegms(z:end, :) = [];%delete excess
setappdata(evi, 'extSegms', extSegms);%external stimulus sent during experiment
    
CutSpnt_Callback([], [], handles) %creat segment and plot

set(evi, 'WindowButtonDownFcn', {@Eview_WindowButtonDownFcn, handles})%reaction on button down
set(evi, 'WindowButtonMotionFcn', {@Eview_WindowButtonMotionFcn, handles})%reaction on button motion
set(evi, 'WindowButtonUpFcn', {@Eview_WindowButtonUpFcn, handles})%reaction on button up




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
manlDet = getappdata(evi, 'manlDet');%manually detected events

if isequal(get(hObject, 'Tag'), 'CutSpnt') %source
    set(handles.CutManual, 'Value', 0);%switch off manually detected events
else
    set(handles.CutSpnt, 'Value', 0);%switch off spontaneous cutting
end
if (~isempty(manlDet) && ~get(handles.CutManual, 'Value')) %try save manlDet
    ynBtn = questdlg('Save manually detected?', 'Save events', 'Yes', 'No', 'No');%ask what to do
    if isequal(ynBtn, 'Yes') %save requested
        pth = getappdata(evi, 'pth');%directory
        flNm = getappdata(evi, 'flNm');%filename

        z = 0;%file number
        while exist([pth, flNm, '_mdet', num2str(z), '.mat'], 'file') %file exit
            z = z + 1;%next file (with other events)
        end
        save([pth, flNm, '_mdet', num2str(z), '.mat'], 'manlDet') %save manually detected events
    end
    setappdata(evi, 'manlDet', []);%manually detected events
    set(handles.CutManual, 'String', 'mdet');%number of manually detected events
end
if (isempty(manlDet) && get(handles.CutManual, 'Value')) %try open manlDet
    ynBtn = questdlg('Load manually detected?', 'Load events', 'Yes', 'No', 'No');%ask what to do
    if isequal(ynBtn, 'Yes') %load requested
        pth = getappdata(evi, 'pth');%directory
        flNm = getappdata(evi, 'flNm');%filename
        [flNm, pth] = uigetfile({'*.mat', 'mat'}, 'Select file', [pth, flNm, '.mat']);%open dialog
        if ~(isequal(flNm, 0) && isequal(pth, 0)) %file was choosed
            varInMatFile = who(matfile([pth, flNm]));%variable in the choosed file
            if ismember('manlDet', varInMatFile)
                load([pth, flNm], 'manlDet');%load manually detected events
                setappdata(evi, 'manlDet', manlDet);%manually detected events
                set(handles.CutManual, 'String', ['mdet ', num2str(length(manlDet))]);
            else
                set(handles.CutManual, 'Value', 0);%switch off manually detected events
            end
        else
            set(handles.CutManual, 'Value', 0);%switch off manually detected events
        end
    else
        set(handles.CutManual, 'Value', 0);%switch off manually detected events
    end
end
extSegms = getappdata(evi, 'extSegms');%external stimulus sent during experiment

%= creat segments of the recordation =%
%segms(:, 1) - position of synchro-event (ms)
%segms(:, 2) - number of channel where syncro-event was detected
%segms(:, 3) - number of sweeps where syncro-event was detected
%segms(:, 4) - number of stimulus (within sweep, in matrix zavp.realStim) or number of trough in brst matrix
%segms(:, 5) - sweep start (seconds from beginning of recording)
if ((get(handles.CutSpnt, 'Value') || isempty(extSegms)) && ~(get(handles.CutManual, 'Value') && ~isempty(manlDet))) %spontwise cutting choosed
    aftStimStr = get(handles.AftStim, 'String');%string value
    if (isequal(aftStimStr, 'NaN') || isempty(aftStimStr)) %not a number
        tW = 10e3;%time window (ms)
    else
        tW = str2double(aftStimStr);%time window (ms)
    end
    if (hd.lActualEpisodes <= 1) %single sweeps (gap-free mode)
        segms = zeros(floor(size(lfp, 1) / tW), 4);%tW-seconds segments
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
            segms(z, 4) = NaN;%number of stimulus (in matrix zavp.realStim) or number of trough in brst matrix
            z = z + 1;
        end
        segms(z:end, :) = [];%delete excess
    else %many sweeps (continuous mode)
        segms = zeros(hd.lActualEpisodes, 4);%segments
        for sw = 1:hd.lActualEpisodes %run over sweeps
            segms(sw, 1) = round(size(lfp, 1) / 2);%position of synchro-event (ms)
            segms(sw, 2) = NaN;%number of channel where syncro-event was detected
            segms(sw, 3) = sw;%number of sweep where syncro-event was detected
            segms(sw, 4) = 1;%number of stimulus (in matrix zavp.realStim) or number of trough in brst matrix
        end
    end
elseif (get(handles.CutManual, 'Value') && ~isempty(manlDet)) %cut by manually detected events
    segms = zeros(length(manlDet), 4);%preallocation of memory for events moments
    for t = 1:length(manlDet) %run over manually detected events
        segms(t, 1) = manlDet(t).t(1);%position of synchro-event (ms from record begin (from 0))
        segms(t, 2) = manlDet(t).ch(1);%number of channel where syncro-event was detected
        segms(t, 3) = manlDet(t).sw;%number of sweep where syncro-event was detected
        segms(t, 4) = t;%number of stimulus (in matrix zavp.realStim) or number of trough in brst matrix
    end
else %cut with respect to external synchronization (created during experiment)
    segms = extSegms;%external stimulus sent during experiment
    tW = min(1e5, round(size(lfp, 1) ./ size(segms, 1)));%length of view window (ms)
end

%= creat time-tags =%
timeStr = repmat(' ', size(segms, 1), 8);%"Moscow" time of the segment begin (part of day)
for t = 1:size(segms, 1) %run over segments
    segms(t, 5) = 0;%sweep start (seconds from beginning of recording)
    if isfield(hd, 'sweepStartInPts') %sweep start exist
        segms(t, 5) = hd.sweepStartInPts(segms(t, 3)) * hd.si / 1e6;%sweep start (seconds from beginning of recording)
    end
    timeStr(t, :) = MoscowTime(hd.recTime(1) + segms(t, 5) + (segms(t, 1) * 1e-3));%"Moscow" time of the segment begin (part of day)
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

lfpLn = getappdata(evi, 'lfpLn');%LFP line handles
spkBr = getappdata(evi, 'spkLn');%spke bar handles
spkDnLn = getappdata(evi, 'spkDnLn');%spike density line handles
clrMapH = getappdata(evi, 'clrMapH');%colormap handles (CSD or MUA)

if isfield(zavp, 'chOrd')
    lfp = lfp(:, zavp.chOrd, :);
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

% cla(lfp_ax)%clear LFP-axis
dltGrphObj = getappdata(evi, 'dltGrphObj');%deletable objects
for t = 1:length(dltGrphObj) %run over object to be deleted
    delete(dltGrphObj(t));%delet object
end
dltGrphObj = [];%no objects to be deleted
setappdata(evi, 'dltGrphObj', []);

for t = 1:length(lfpLn) %run over objects
    set(lfpLn(t), 'Visible', 'off')%inivisible now
end
for t = 1:length(spkBr) %run over objects
    set(spkBr(t), 'Visible', 'off')%inivisible now
end
for t = 1:length(spkDnLn) %run over objects
    set(spkDnLn(t), 'Visible', 'off')%inivisible now
end
for t = 1:length(clrMapH) %run over objects
    set(clrMapH(t), 'Visible', 'off')%inivisible now
end

setappdata(evi, 'dFrmH', []);%drag frame handle
setappdata(evi, 'bgnXY', []);%begin point of rectangle draw
setappdata(evi, 'endXY', []);%end point of rectangle draw
chnlStr = cell(0);%vertical axis labels (channels)

%= run over channel groups and plot traces individual settings =%
unqChnlGrp = unique(horzcat(chTabData{:, 4}));%names of unique groups
initLvl = 0;%initial drawing level
for z = 1:length(chanSettingGrp) %run over all channels
    chanSettingGrp(z).realsernum = 0;%actual serial number of the channel on a graph
end
realsernum = 1;%actual serial number of the first visible channel
for cg = unqChnlGrp %run over unique groups
    jj = horzcat(chanSettingGrp(:).grpName);%all groups
    cgch = find(jj == cg);%current group channels
    rCh = find(horzcat(chTabData{:, 3}));%all selected channels
    rCh = intersect(rCh, cgch);%number of visible (selected) channels of current group rCh array
    for z = 1:length(rCh) %run over selected channels of the group
        chanSettingGrp(rCh(z)).realsernum = realsernum;%actual serial number of the channel on a graph
        realsernum = realsernum + 1;%actual serial number of next visible channel
    end
    
    if ~isempty(rCh) %any channels was selected
    rawData = chanSettingGrp(rCh(1)).rawSignal;%load raw data if true
    
    %= timevector for raw and resampled signals =%
    if rawData %raw data requested
        prm.Fs = 1 / zavp.siS;%sampling frequency for raw data
    else %resampled data requested
        prm.Fs = zavp.dwnSmplFrq;%sampling frequency for rare data
    end
    tm = segmEdge(1):(1e3 / prm.Fs):segmEdge(2);%time matrix for raw or downsampled data (ms)
    if ((segmEdge(2) - tm(end)) >= ((1e3 / prm.Fs) / 2))
        tm = [tm, tm(end) + (1e3 / prm.Fs)];
    end
    %= end of (timevector for raw and resampled signals) =%
    

    %= common reference (step 1) =%
    if chanSettingGrp(rCh(1)).comRef %common reference subtraction
        tmpChnl = cgch;%requested channels
    else %without common reference subtraction
        tmpChnl = rCh;%requested channels
    end
    if rawData %raw data is requested
        fLFP = ZavSynchLFP(zavp, hd, segms(sn, :), segmEdge, lfp, tmpChnl, rawData, nlxVer);%lfp phased with respect to stimulus moments
    else %resampled data is requested
        fLFP = lfp(:, tmpChnl, :);%no filters
    end
    %= end of (common reference (step 1)) =%
    
        
    %= LFP-segmentation =%
    if rawData %raw data is requested
        whtLFP = fLFP;%lfp phased with respect to stimulus moments
    else %resampled data is requested
        whtLFP = ZavSynchLFP(zavp, hd, segms(sn, :), segmEdge, fLFP, 1:size(fLFP, 2), rawData, nlxVer);%lfp phased with respect to stimulus moments
    end
    %= end of (LFP-segmentation) =%
    

    %= DC-shift compensation =%
    if chanSettingGrp(rCh(1)).dcShift %DC compensation
        ii = 1:min(20, size(whtLFP, 1));%time near segment begin
        %ii = (tm > -1e1) & (tm < 1e1);%time near synchro-event
        ajj = median(whtLFP(ii, :, :), 1);%base-level
        %ajj = median(whtLFP(:, :, :), 1);%base-level
        whtLFP = whtLFP - repmat(ajj, [size(whtLFP, 1), 1, 1]);%subtract DC component
    end
    %= end of (DC-shift compensation) =%

    %= common reference (step 2) =%
    if chanSettingGrp(rCh(1)).comRef %common reference subtraction
        ii = 1:min(20, size(whtLFP, 1));%time near segment begin
        %ii = (tm > -1e1) & (tm < 1e1);%time near synchro-event
        ajj = median(whtLFP(ii, :, :), 1);%base-level
        wLFPcomm = mean(whtLFP - repmat(ajj, [size(whtLFP, 1), 1, 1]), 2);%common reference with adduction to zero-level of each channel in a group
        [~, irch] = intersect(cgch, rCh);%number of selected (visible) channel of the group
        whtLFP = whtLFP(:, irch, :) - repmat(wLFPcomm, [1, length(rCh), 1]);%subtract common reference
    end
    allSwLFP = whtLFP;%all sweeps lfp
    whtLFP = mean(whtLFP, 3);%average by segments (after transformations)
    %= end of (common reference (step 2)) =%
    

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
    %= end of (time axis parameters) =%
    
    
    %= drawing step (plot acceleration) =%
    rarStep = 1;% + floor(size(whtLFP, 1) / 10e3);%plot with step
    if (rarStep > 1)
        ajj = whtLFP;%copy of LFP segment
        whtLFP = zeros(length(resample(ajj(:, 1), 1, rarStep)), size(ajj, 2));
        for t = 1:size(whtLFP, 2) %run over channels
            whtLFP(:, t) = resample(ajj(:, t), 1, rarStep);%downsampled LFP
        end
        ajj = allSwLFP;%copy of LFP segment
        allSwLFP = zeros(length(resample(ajj(:, 1, 1), 1, rarStep)), size(ajj, 2), size(ajj, 3));
        for z = 1:size(ajj, 3) %run over segments
            for t = 1:size(allSwLFP, 2) %run over channels
                allSwLFP(:, t, z) = resample(ajj(:, t, z), 1, rarStep);%downsampled LFP
            end
        end
        tm = interp1(1:length(tm), tm, linspace(1, length(tm), size(whtLFP, 1)));%downsampled time
    end
    %= end of (drawing step (plot acceleration)) =%
    

    %= plot colormapped CSD =%
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
        
        %spatial filtration
        lfpCSD = filtfilt([0.607, 1, 0.607], 2.213, lfpCSD')';%3-channel spatial smoothing as in elephy
%         for t = 1:size(lfpCSD, 1) %run over time
%             lfpCSD(t, :) = smooth(lfpCSD(t, :), 5);%spatial filtration
%         end

        %interpolation (along channels)
        ajj = lfpCSD;%copy of original CSD map
        lfpCSD = zeros(size(ajj, 1), length(1:0.2:size(lfpCSD, 2)));
        for t = 1:size(lfpCSD, 1) %run over time
            lfpCSD(t, :) = interp1(1:(size(whtLFP, 2) - 2), ajj(t, :), 1:0.2:(size(whtLFP, 2) - 2));%interpolation along channels
        end

        %subplot(lfp_ax);%set lfp_ax as a parent for next plot
        %imagesc(tm / timeScalCoeff, initLvl + (2:(length(rCh) - 1)), lfpCSD');
        set(clrMapH(cg), 'XData', tm / timeScalCoeff, 'YData', initLvl + (2:(length(rCh) - 1)), 'CData', lfpCSD', 'Visible', 'on')
        caxis(lfp_ax, [str2double(get(handles.CSDcMapMin, 'String')), str2double(get(handles.CSDcMapMax, 'String'))])%colorbar range
        
        colormap('jet')
    end
    %= end of (plot CSD) =%
    

    %= plot spikes density =%
    if (chanSettingGrp(rCh(1)).pltSpkMap || chanSettingGrp(rCh(1)).pltMUAfreq) %plot map of spikes frequency
        binStep = [chanSettingGrp(rCh(1)).muaBin, chanSettingGrp(rCh(1)).precisF];%[bin size (ms), precision (must be integer)]
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

        dCoeff(1) = -1; dCoeff(2) = chanSettingGrp(rCh(1)).muaScale;%characteristic sizes
        dCoeff(3) = dCoeff(1) * dCoeff(2);%scaling coefficient

        clr = [0, 0, 1];%color of spikes density line
        if chanSettingGrp(rCh(1)).pltSpkMap %plot colormapped MUA density
        %     %interpolation (smooth along channels)
        %     spkDensSC = zeros(size(spkDens, 1), length(1:0.25:length(rCh)));
        %     for t = 1:size(spkDens, 1) %run over time
        %         spkDensSC(t, :) = interp1(1:length(rCh), spkDens(t, :), 1:0.25:length(rCh));%interpolation along channels
        %     end
            spkDensSC = spkDens;
            %smooth along time
            for ch = 1:size(spkDens, 2) %rum over channels
                spkDensSC(:, ch) = smooth(spkDens(:, ch), 3);%smooth along time
            end

            %suplot(lfp_ax);%set lfp_ax as a parent for next plot
            %imagesc(basbins / timeScalCoeff, initLvl + (1:length(rCh)), spkDensSC');
            set(clrMapH(cg), 'XData', basbins / timeScalCoeff, 'YData', initLvl + (1:length(rCh)), 'CData', spkDensSC', 'Visible', 'on')
            shading flat;
            caxis(lfp_ax, [str2double(get(handles.SpkCMapMin, 'String')), str2double(get(handles.SpkCMapMax, 'String'))])%colorbar range
            clr = [1, 1, 1];%color of spikes density line
            
            colormap('jet')
        end %if ~chanSettingGrp(rCh(1)).pltCSDmap %no CSD plot choosed

        %= plot spikes densities
        if chanSettingGrp(rCh(1)).pltMUAfreq %plot MUA frequency traces
            %linHandl = plot(lfp_ax, basbins / timeScalCoeff, (spkDens / dCoeff(3)) + repmat(initLvl + (1:length(rCh)), length(basbins), 1), clr, 'LineWidth', 1);
            %set(linHandl, 'ButtonDownFcn', {@Line_ButtonDownFcn, handles})%reaction on mouse-click
            for ch = 1:size(spkDens, 2) %run over channels
                set(spkDnLn(initLvl + ch), 'XData', basbins / timeScalCoeff, 'YData', (spkDens(:, ch) / dCoeff(3)) + initLvl + ch, 'Color', clr, 'Visible', 'on')
                set(spkDnLn(initLvl + ch), 'UserData', {dCoeff, ch, timeScalCoeff, initLvl, rCh(ch)});%save amplitude coefficient and offset
            end

    %         [~, p1] = ds2nfu([0, 0], [0, abs(1 / dCoeff(1))]);%diff(p1) = length of scale bar in normalized units
    %         set(handles.MUAscb, 'Position', [0.525, 0.15 0, diff(p1)], 'LineWidth', 5, 'Units', 'normalized', 'TextRotation', 90, 'HeadStyle', 'none',...
    %             'String', [num2str(round(dCoeff(2) * 10) / 10), 'spikes/', num2str(binStep(1)), 'ms'], 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');
        end
    end
    %= end of (plot spikes density) =%
    
    
    %= plot LFP =%
    if chanSettingGrp(rCh(1)).pltLFP %plot LFP traces
        %= LFP scaling =%
        dCoeff(1) = -1.0; dCoeff(2) = chanSettingGrp(rCh(1)).lfpScale;%characteristic sizes (amplitude coefficient)
        dCoeff(3) = dCoeff(1) * dCoeff(2) * 1e3;%scaling coefficient
        whtLFP = whtLFP / dCoeff(3);%scaled signals
        allSwLFP = allSwLFP / dCoeff(3);%all sweeps lfp
    
        if (size(allSwLFP, 3) > 1) %averaging
            for ch = 1:size(whtLFP, 2) %run over channels
                for z = 1:size(allSwLFP, 3) %run over segments
                    h = plot(lfp_ax, tm / timeScalCoeff, allSwLFP(:, ch, z) + initLvl + ch, 'Color', 0.9 * ones(1, 3));
                    dltGrphObj = [dltGrphObj; h];%add object to be deleted
                end
            end
        end
        
        %linHandl = plot(lfp_ax, tm / timeScalCoeff, whtLFP + repmat(initLvl + (1:length(rCh)), length(tm), 1), 'k', 'LineWidth', 1);
        %set(linHandl, 'ButtonDownFcn', {@Line_ButtonDownFcn, handles})%reaction on mouse-click
        for ch = 1:size(whtLFP, 2) %run over channels
            set(lfpLn(initLvl + ch), 'XData', tm / timeScalCoeff, 'YData', whtLFP(:, ch) + initLvl + ch, 'Visible', 'on', 'Color', zeros(1, 3))
            set(lfpLn(initLvl + ch), 'UserData', {dCoeff, ch, timeScalCoeff, initLvl, rCh(ch)});%save amplitude coefficient and offset
        end

        %= plot bursts =%
        if (~isempty(brst) && (length(sn) == 1) && isequal(get(handles.ShowBrst, 'State'), 'on')) %single segment
            jj = vertcat(brst(:).t);
            p1 = find((jj(:, 1) < (segms(sn, 1) + segmEdge(2))) & (jj(:, 2) > (segms(sn, 1) + segmEdge(1))));%find visible bursts
            for t = p1' %run over visible bursts
                [~, ii] = intersect(rCh, brst(t).ch);%common channels
                ii = ii(:)';
                if ~isempty(ii) %some channels contain the burst
                    jj = ((tm >= (brst(t).t(1) - segms(sn, 1))) & (tm <= (brst(t).t(2) - segms(sn, 1))));
                    h = plot(lfp_ax, tm(jj) / timeScalCoeff, whtLFP(jj, ii) + repmat(initLvl + ii, sum(jj), 1), 'm', 'LineWidth', 1);
                    dltGrphObj = [dltGrphObj; h];%add object to be deleted
                end
            end
        end
    end
    %= end of (plot LFP) =%
    

    %= plot spike bars =%
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
            %plot(lfp_ax, spkCdens(:, 1) / timeScalCoeff, spkCdens(:, 2), 'r', 'LineWidth', 1.5)
            set(spkBr(initLvl + ch), 'XData', spkCdens(:, 1) / timeScalCoeff, 'YData', spkCdens(:, 2), 'Visible', 'on')
        end
    end
    %= end of (plot spike bars) =%
    

    % %= scale bars =% not ended
    % p1 = str2double(get(handles.LFPScale, 'String'));
    % if ~isempty(p1)
    %     if ~isfinite(p1)
    %         p1 = max(max(abs(whtLFP)));
    %         set(handles.LFPScale, 'String', num2str(p1))
    %     end
    % end
    %
    % set(get(lfp_ax, 'YLabel'), 'String', 'channels', 'Units', 'normalized', 'Position', [-0.06, 0.6, 1])
    % [~, p1] = ds2nfu([0, 0], [0, abs(1 / dCoeff(1))]);%diff(p1) = length of scale bar in normalized units
    % set(handles.LFPscb, 'Position', [0.025, 0.02, 0, diff(p1)], 'LineWidth', 5, 'Units', 'normalized', 'TextRotation', 90, 'HeadStyle', 'none',...
    %     'String', [num2str(round(dCoeff(2))), ' \muV'], 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', 'FontSize', 9);
    %
    % %= end of (scale bars) =%
    
    %= vertical axis labels =%
    for t = 1:length(rCh)
        if ~isempty(chTabData{rCh(t), 2}) %alternames was signed
            chnlStr = [chnlStr; chTabData(rCh(t), 2)];%alter name of channels
        else
            chnlStr = [chnlStr; chTabData(rCh(t), 1)];%numbered channels
        end
    end
    %= end of (vertical axis labels) =%

    initLvl = initLvl + length(rCh);%drawing level for next group
    end %end of (if ~isempty(rCh) %any channels was selected)
end %end of (for cg = unqChnlGrp %run over unique groups)
setappdata(evi, 'chanSettingGrp', chanSettingGrp);%save channels group setting

if (size(allSwLFP, 3) > 1) %average requested
    h = get(lfp_ax, 'Children');
    jj = zeros(length(h), 1);
    for t = 1:length(h) %run over graphic objects
        if isequal(get(h(t), 'Tag'), 'lfpLn') %if lfp-line object detected (use findobj?)
            jj(t) = 1;
        elseif isequal(get(h(t), 'Tag'), 'clrMap') %if color map object detected (use findobj?)
            jj(t) = 2;
        end
    end
    set(lfp_ax, 'Children', [h(jj == 1); h(jj == 0); h(jj == 2)]);%set color maps to background
end

%= complete graph creation =%
rCh = find(horzcat(chTabData{:, 3}));%all selected channels

y_lmt = [-0.05, length(rCh) + 1.05];%vertical limits
set(lfp_ax, 'YDir', 'reverse', 'XLim', segmEdge / timeScalCoeff, 'YLim', y_lmt)
set(get(lfp_ax, 'XLabel'), 'String', ['time, ', timeUnit])
set(lfp_ax, 'YTick', 1:length(chnlStr), 'YTickLabel', chnlStr)
if get(handles.ShowMarkers, 'Value') %show marker
    h = plot(lfp_ax, [0, 0], y_lmt, 'c', 'LineWidth', 2);%nul-time line (synchro-event)
    dltGrphObj = [dltGrphObj; h];%add object to be deleted
end

if (length(sn) == 1) %single segment is shown
    set(evi, 'Name', ['Eview - ', flNm, ' (', timeStr(sn, :), ')'])
    h = text(0, y_lmt(1) + diff(y_lmt) / 80, timeStr(sn, :));
    dltGrphObj = [dltGrphObj; h];%add object to be deleted
    
    if get(handles.ShowMarkers, 'Value') %show marker
        extSegms = getappdata(evi, 'extSegms');%external stimulus sent during experiment
        for t = 1:size(extSegms, 1) %run over external stimulus
            if ((extSegms(t, 1) > (segms(sn, 1) + segmEdge(1))) && (extSegms(t, 1) < (segms(sn, 1) + segmEdge(2))) && (extSegms(t, 3) == segms(sn, 3))) %the event is in visible range
                h = plot(lfp_ax, (extSegms(t, 1) - extSegms(sn, 1)) * [1, 1] / timeScalCoeff, y_lmt, 'c', 'LineWidth', 2);%event-line
                dltGrphObj = [dltGrphObj; h];%add object to be deleted
            end
        end

        for t = 2:length(eventsTags) %run over events
            if ((eventsTags{t, 1} > (segms(sn, 1) + segmEdge(1))) && (eventsTags{t, 1} < (segms(sn, 1) + segmEdge(2)))) %the event is in visible range
                h = plot(lfp_ax, (eventsTags{t, 1} - segms(sn, 1)) * [1, 1] / timeScalCoeff, y_lmt, 'g', 'LineWidth', 2);%event-line
                dltGrphObj = [dltGrphObj; h];%add object to be deleted
                h = text((eventsTags{t, 1} - segms(sn, 1)) / timeScalCoeff, y_lmt(1) + diff(y_lmt) / 80, eventsTags{t, 2});
                dltGrphObj = [dltGrphObj; h];%add object to be deleted
            end
        end
    end
    
    manlDet = getappdata(evi, 'manlDet');%manually detected events
    for t = 1:length(manlDet) %run over manually detected events
        if ((manlDet(t).t > (segms(sn, 1) + segmEdge(1))) && (manlDet(t).t < (segms(sn, 1) + segmEdge(2))) && (manlDet(t).sw == segms(sn, 3))) %the event is in visible range
            h = plot(lfp_ax, (manlDet(t).t - segms(sn, 1)) * [1, 1] / timeScalCoeff, y_lmt, 'b', 'LineWidth', 1.5);%event-line
            dltGrphObj = [dltGrphObj; h];%add object to be deleted
            h = text((manlDet(t).t - segms(sn, 1)) / timeScalCoeff, y_lmt(1) + diff(y_lmt) / 80, 'md');
            dltGrphObj = [dltGrphObj; h];%add object to be deleted
        end
    end
else
    set(evi, 'Name', ['Eview - ', flNm, ' (', timeStr(sn(1), :), ' - ', timeStr(sn(end), :), ')'])
end
%= end of (complete graph creation) =%
setappdata(evi, 'dltGrphObj', dltGrphObj)%new array with deletable object handles



% --- Executes on button press in Average.
function Average_Callback(hObject, eventdata, handles)
% hObject    handle to Average (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(hObject, 'Value', 1);
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
    'String', {['X:', num2str(round(xy(1) * 10) / 10)], ['Y:', num2str(round(xy(2) * 100) / 100)], timeStr(m, :)}, ...
    'BackgroundColor', [1, 1, 1], 'EdgeColor', [0, 1, 0], 'UserData', {xy, hObject}, 'FontSize', 8);
dataCurs(nCrs).l = annotation('line', 'Units', 'normalized', 'Position', [p1, 1 - p2 - 0.02, 0.0, -0.02], ...
    'Color', [0, 1, 0]);

if (nCrs > 1) %many labels
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

lfpLn = getappdata(evi, 'lfpLn');%lfp-line handles
if ~isempty(lfpLn) %data was loaded
    
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
        TunChannels(floor(bgnXY(1, 2)), handles)%call channels tuning dialog
    end
    return;%exit from Eview_WindowButtonDownFcn
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
end %if ~isempty(lfpLn) %data was loaded


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
axFrm = getappdata(evi, 'axFrm');%handle of axes with frame
endXY = get(axFrm, 'CurrentPoint');%right lower point
if (~isempty(bgnXY) && ~isempty(endXY)) %click-data exist
    x_lim = sort([bgnXY(1, 1), endXY(1, 1)]);%horizontal range
    y_lim = sort([bgnXY(1, 2), endXY(1, 2)]);%vertical range
else
    return;%exit from Eview_WindowButtonUpFcn
end

%immediately remove the frame
delete(getappdata(evi, 'dFrmH'))
setappdata(evi, 'dFrmH', []);%drag frame handle
setappdata(evi, 'bgnXY', []);%begin point of rectangle draw
setappdata(evi, 'endXY', []);%end point of rectangle draw

hd = getappdata(evi, 'hd');%header
segms = getappdata(evi, 'segms');%segments
sn = str2double(get(handles.CurrSegm, 'String'));%segments of signal (sweeps)

lfpLn = getappdata(evi, 'lfpLn');%lfp-line handles
if ~isempty(lfpLn) %data was loaded
    dCoeff = get(lfpLn(1), 'UserData');%ampl-time-scale factors

    %selected channels
    chanSettingGrp = getappdata(evi, 'chanSettingGrp');%load channels group setting
    rCh = 1:hd.nADCNumChannels;%all possible channels
    visChNum = horzcat(chanSettingGrp(:).realsernum);%actual serial numbers of the channel on a graph
    rCh = rCh(visChNum > 0);%visible channels only
    visChNum(visChNum == 0) = [];%delete zeros
    ii = rCh;%copy
    for z = 1:length(visChNum) %run over visible channels
       rCh(visChNum(z)) = ii(z);%actual order of channels
    end

    if ((diff(x_lim) > 0) && (diff(y_lim) > 0)) %drag was done
        bT.t = segms(sn, 1) + x_lim * dCoeff{3};%begin-end of currently defined burst (or highlighted area)
        jj = ceil(y_lim(1)):floor(y_lim(2));%highlighted vertical range
        jj = jj((jj > 0) & (jj <= length(rCh)));%selected channels
        bT.ch = rCh(jj);%selected channels
        bT.t(1) = max([1, bT.t(1)]);%left edge of the frame
        bT.t(2) = min([bT.t(2), (diff(hd.recTime) * 1e3) - 1]);%right edge of the frame
        if isequal(get(handles.PwrSpectrWin, 'State'), 'on') %spectrum in window selected
            %plot spectra
            if (~isempty(bT.ch) && (diff(bT.t) > 5)) %traces selected
                if (get(handles.WinAvr, 'Value')) %average segments
                    z = str2double(get(handles.CurrSegm, 'String'));%first segment to be averaged
                    t = str2double(get(handles.AvrLen, 'String'));%averaging length (number of segments to be averaged)
                    sn = max(1, z):min((z + t - 1), size(segms, 1));%selected segments
                    %sn = 1:size(segms, 1);%all segments
                end
                %if (chanSettingGrp(slctCh).pltLFP) %LFP is shown
                    lfpSpctFigH = figure;
                    set(lfpSpctFigH, 'Name', 'LFP spectra', 'NumberTitle', 'off', 'Units', 'normalized');
                    lfp = getappdata(evi, 'lfp');
                    zavp = getappdata(evi, 'zavp');
                    hd = getappdata(evi, 'hd');
                    whtLFP = ZavSynchLFP(zavp, hd, [segms(sn, 1) + x_lim(1) * dCoeff{3}, repmat([NaN, 1, 1], numel(sn), 1)], [0, round(diff(bT.t))], lfp, bT.ch, 0, 0);
                    whtLFP = ZavFilter(whtLFP, 1e3, 'high', 5, 2);
                    whtLFP = mean(whtLFP, 3);

                    prm.Fs = zavp.dwnSmplFrq;%discretization frequency for signals to be treated (Hz)
                    prm.fpass = [1, 100];%frequency range (predominantly for sensory response)
                    prm.tapers = [diff(prm.fpass) * round(diff(bT.t)) * 1e-3, 5];%[3, 5];%tapers
                    prm.pad = 2;%padding factor
                    prm.trialave = 0;%channel by channel, sweep by sweep (don't average when 0)

                    [pwrSpct, frq] = mtspectrumc(whtLFP, prm);%mean power spectrum (mean by sweeps)
                    hold on
                    plot(frq, pwrSpct)
                    plot(frq, mean(pwrSpct, 2), 'k', 'LineWidth', 2)
                    xlim(frq([1, end]))
                    xlabel('frequency, Hz')
                    ylabel('power')
                    lgndStr = {};
                    for z = bT.ch
                        lgndStr{end + 1} = ['ch', num2str(z)];
                    end
                    lgndStr{end + 1} = 'average';
                    legend(lgndStr)
                    set(lfpSpctFigH, 'UserData', bT.t / 1e3)%seconds from record begin
                    expBtnH = uicontrol('Parent', lfpSpctFigH, 'Units', 'normalized', 'Style', 'pushbutton', 'String', 'export', ...
                                        'Position', [0.82, 0.935, 0.1, 0.05], ...
                                        'Callback', ['h = get(gca, ''Children''); frq = get(h(1), ''XData'')''; pwrSpct = []; ', ...
                                                     'for t = 2:(length(h) - 0), pwrSpct = [pwrSpct, get(h(t), ''YData'')'']; end; ', ...
                                                     'tRang = get(gcf, ''UserData''); ', ...
                                                     'disp(''export done'')'], ... 'save(''eviewSpectra'', ''frq'', ''pwrSpct'', ''tRang'')'
                                        'BackgroundColor', [0.831, 0.816, 0.784]);
                %end

                if ~isempty(spks) %spikes exist any(chanSettingGrp(bT.ch).pltMUAfreq) %MUA trace is shown
                    muaSpctFigH = figure;
                    set(muaSpctFigH, 'Name', 'MUA spectra', 'NumberTitle', 'off', 'Units', 'normalized')
                    spks = getappdata(evi, 'spks');
                    binStep = [str2double(get(handles.SpkBin, 'String')), str2double(get(handles.HistPrecision, 'String'))];%[bin size (ms), precision (must be integer)]
                    [~, basbins] = ZavOverHist(spks(bT.ch(1), 1).tStamp, binStep, bT.t);%fast overlaped histogram (histc)
                    spkDens = zeros(length(basbins), numel(bT.ch), numel(sn));%memory preallocation for density of spikes in time (histogram)

                    for z = 1:numel(sn) %run over requested segments
                        jj = segms(sn(z), 1) + x_lim * dCoeff{3};%time segment boundary (ms from begin of file)
                        sw = segms(sn(z), 3);%sweep
                        for ch = 1:numel(bT.ch) %run over channels
                            spkDens(:, ch, z) = ZavOverHist(spks(bT.ch(ch), sw).tStamp, binStep, jj);%fast overlaped histogram (histc)
                        end
                    end
                    spkDens = mean(spkDens, 3) / binStep(1);%mean spikes density (spikes / ms)

                    prm.Fs = binStep(2) / (binStep(1) * 1e-3);%discretization frequency for signals to be treated (Hz)
                    prm.fpass = [10, min(prm.Fs / 2, 1000)];%frequency range (predominantly for sensory response)
                    prm.tapers = [diff(prm.fpass) * round(diff(bT.t)) * 1e-3, 5];%[3, 5];%tapers
                    prm.pad = 2;%padding factor
                    prm.trialave = 0;%channel by channel, sweep by sweep (don't average when 0)

                    spkDens = ZavFilter(spkDens, prm.Fs, 'high', 50, 2);%high cut filtration

                    [pwrSpct, frq] = mtspectrumc(spkDens, prm);%mean power spectrum (mean by sweeps)
                    hold on
                    plot(frq, pwrSpct)
                    plot(frq, mean(pwrSpct, 2), 'k', 'LineWidth', 2)
                    xlim(frq([1, end]))
                    xlabel('frequency, Hz')
                    ylabel('power')
                    lgndStr = {};
                    for z = bT.ch
                        lgndStr{end + 1} = ['ch', num2str(z)];
                    end
                    lgndStr{end + 1} = 'average';
                    legend(lgndStr)
                    set(muaSpctFigH, 'UserData', bT.t / 1e3)%seconds from record begin
                    expBtnH = uicontrol('Parent', muaSpctFigH, 'Units', 'normalized', 'Style', 'pushbutton', 'String', 'export', ...
                                        'Position', [0.88, 0.92, 0.1, 0.05], ...
                                        'Callback', ['h = get(gca, ''Children''); frqSpk = get(h(1), ''XData'')''; pwrSpctSpk = []; ', ...
                                                     'for t = 2:(length(h) - 0), pwrSpctSpk = [pwrSpctSpk, get(h(t), ''YData'')'']; end; ', ...
                                                     'tRang = get(gcf, ''UserData''); ', ...
                                                     'disp(''export done'')'], ... 'save(''eviewSpectra'', ''frq'', ''pwrSpct'', ''tRang'')'
                                        'BackgroundColor', [0.831, 0.816, 0.784]);
                end %end of (if ~isempty(spks) %spikes exist)
            end %end of (if (~isempty(bT.ch) && (diff(bT.t) > 5)) %traces selected)
        else %complete selected burst
            brst = getappdata(evi, 'brst');
            if ~isempty(brst) %bursts exist
                z = length(brst) + 1;%number of the new burst
            else %no one burst yet
                clear brst
                brst = struct('t', [], 'ch', []);%creat bursts data
                z = 1;%number of the new burst
            end

            sb = false;%if any burst exist on selected area
            b = [];%number of bursts within selected area
            jj = vertcat(brst(:).t);%all detected bursts
            if ~isempty(jj) %bursts exist
                b = find((jj(:, 1) <= bT.t(1)) & (jj(:, 2) >= bT.t(1)));%left board on burst
                b = [b; find((jj(:, 1) <= bT.t(2)) & (jj(:, 2) >= bT.t(2)))];%right board on burst
                b = [b; find((jj(:, 1) >= bT.t(1)) & (jj(:, 1) <= bT.t(2)))];%number of bursts within selected area
                b = [b; find((jj(:, 2) >= bT.t(1)) & (jj(:, 2) <= bT.t(2)))];%number of bursts within selected area
                b = unique(b);%unique bursts numbers
                ii = [];
                for t = 1:length(b)
                    if ~isempty(intersect(brst(b(t)).ch, bT.ch)) %common channels exist
                        ii = [ii, t];
                    end
                end
                b = b(ii);%take channels into account
                if ~isempty(b) %burst exist within selected area
                    sb = true;%some burst exist on selected range
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
    elseif (isequal(get(handles.ManualDetect, 'State'), 'on') || ... manual detection choosed
            get(handles.CutManual, 'Value') ... search inside manually detected
           )
        manlDet = getappdata(evi, 'manlDet');%manually detected events
        lfp_ax = handles.LFP_ax;%LFP axis handle
        if isempty(manlDet) %no one manually detected events
            manlDet = struct('t', [], 'ch', [], 'subT', [], 'subCh', [], 'sw', []);%initiate structure for events
            k = 1;%number of manually detected events
        else
            k = length(manlDet) + 1;%number of manually detected events
        end
        
        mdPrm = get(handles.ManDetPrm, 'Data');%manual detection parameters
        srChh = [];%number of needed LFP-handle (main channel)
        subSrChh = [];%number of needed LFP-handle (subchannel)
        for z = 1:length(lfpLn) %run over lines
            if isequal(get(lfpLn(z), 'Visible'), 'on') %channel is drawn
                usrDat = get(lfpLn(z), 'UserData');%get attributes
                if (usrDat{5} == mdPrm{1, 2}) %wanted channel for detection
                    srChh = z;%number of needed LFP-handle
                end
                if (usrDat{5} == mdPrm{4, 2}) %wanted subchannel
                    subSrChh = z;%number of needed LFP-handle
                end
            end
        end
        
        if get(handles.CutManual, 'Value') %search inside manually detected
            srChh = subSrChh;%replace by subchannel
            k = str2double(get(handles.CurrSegm, 'String'));%shown event number
        end
        if (isempty(srChh))
            errordlg('Selected channel is not shown', 'Manual detection error', 'modal')
            return
        end
        
        manlDet(k).ch = mdPrm{1, 2};%main channels for manual detection
        manlDet(k).subCh = mdPrm{4, 2};%subchannel
        
        tmLFP = get(lfpLn(srChh), 'XData') * usrDat{3};%time vector (ms)
        whtLFP = (get(lfpLn(srChh), 'YData') - (usrDat{2} + usrDat{4})) * usrDat{1}(1);%trace
        whtLFP = (2 * (mdPrm{2, 2} > 0) - 1) * whtLFP;%inverse if need
        p1 = find(tmLFP >= round(x_lim(1) * usrDat{3}), 1, 'first');%clicked point
        if (mdPrm{3, 2} > 0) %search range
            tm = -mdPrm{3, 2}:mdPrm{3, 2};%searching range
            tm(((p1 + tm) < 1) | ((p1 + tm) > length(whtLFP))) = [];%exclude outlyers
            p2 = ZavFindMins(whtLFP(p1 + tm));%local minima
            if ~isempty(p2) %extremum found
                [~, t] = min(whtLFP(p1 + tm(p2)));%extremum
                p3 = p1 + tm(p2(t));%point of local minimum (ms from visible segment begin)
            else %no extremum
                errordlg('No extremum. Click other point or change polarity', 'Manual detection error', 'modal')
                return
            end
        else %exact point
            p3 = p1;%exactly clicked point
        end
        
        if ~get(handles.CutManual, 'Value') %manually detected on main channel
            manlDet(k).t = segms(sn, 1) + tmLFP(p3);%absolute time of manually detected event (ms from record start)
            
            %plot manDet-marker
            dltGrphObj = getappdata(evi, 'dltGrphObj');%array with deletable object handles
            y_lmt = [-0.05, length(rCh) + 1.05];%vertical limits
            h = plot(lfp_ax, (manlDet(k).t - segms(sn, 1)) * [1, 1] / usrDat{3}, y_lmt, 'b', 'LineWidth', 1.5);%event-line
            dltGrphObj = [dltGrphObj; h];%add object to be deleted
            h = text((manlDet(k).t - segms(sn, 1)) / usrDat{3}, y_lmt(1) + diff(y_lmt) / 80, 'md');
            dltGrphObj = [dltGrphObj; h];%add object to be deleted
            setappdata(evi, 'dltGrphObj', dltGrphObj)%new array with deletable object handles
        else %search inside manually detected
            manlDet(k).subT = segms(sn, 1) + tmLFP(p3);%absolute time of manually detected event (ms from record start)
            
            %plot subBanDet-marker
            dltGrphObj = getappdata(evi, 'dltGrphObj');%array with deletable object handles
            y_lmt = [-0.05, length(rCh) + 1.05];%vertical limits
            h = plot(lfp_ax, (manlDet(k).subT - segms(sn, 1)) * [1, 1] / usrDat{3}, y_lmt, 'y', 'LineWidth', 1.5);%event-line
            dltGrphObj = [dltGrphObj; h];%add object to be deleted
            h = text((manlDet(k).subT - segms(sn, 1)) / usrDat{3}, y_lmt(1) + diff(y_lmt) / 80, 'md');
            dltGrphObj = [dltGrphObj; h];%add object to be deleted
            setappdata(evi, 'dltGrphObj', dltGrphObj)%new array with deletable object handles
        end
        manlDet(k).sw = segms(sn, 3);%sweep number
        
        %= delete coinciding events =%
        %check time, channel and sweep
        for t = [1:(k - 1), (k + 1):(length(manlDet) - 1)] %run over previous events
            if ((manlDet(t).t == manlDet(k).t) && ...
                (manlDet(t).ch == manlDet(k).ch) && ...
                (manlDet(t).sw == manlDet(k).sw) ...
               )
                manlDet(k) = [];%coincide with previous
                break;
            end
        end
        %= end of (delete coinciding events) =%
        
        jj = vertcat(manlDet(:).t);%time of manually detected events
        [~, ii] = sort(jj(:, 1));%sorted events
        manlDet = manlDet(ii);%sorted events
        
        setappdata(evi, 'manlDet', manlDet);%manually detected events
        set(handles.CutManual, 'String', ['mdet ', num2str(length(manlDet))]);
    end
end %if ~isempty(lfpLn) %data was loaded

        
function TunChannels(slVisCh, handles)
%call dialog for channels tuning
%
%slVisCh - number of selected visible channel
%handles - structure with handles and user data (see GUIDATA)
%
evi = handles.Eview;%handle of the application window
chanSettingGrp = getappdata(evi, 'chanSettingGrp');%load channels group setting

if (slVisCh <= 0) %non integer
    slVisCh = 1;
end
visChNum = horzcat(chanSettingGrp(:).realsernum);%actual serial numbers of the channel on a graph
slctCh = find(visChNum == slVisCh);%real number of channel

%= creat dialog window =%
tunFig = figure('Name', 'Tune channels view', 'Units', 'normalized', 'MenuBar','none', ...
    'NumberTitle', 'off', 'Position', [0.40546875 0.5419921875 0.36796875 0.2392578125], ...
    'Resize', 'off', 'Tag', 'TunChanFig', 'Visible', 'off', 'WindowStyle', 'modal', 'Color', [0.831, 0.816, 0.784]);

tfh(1:36) = struct('h', []);%structure with tune figure control handles

tfh(1).h = uicontrol('Parent', tunFig, 'Units', 'normalized', ...
                    'Style', 'pushbutton', ...
                    'Tag', 'ApplySets', ...
                    'Value', [], ...
                    'String', 'OK', ...
                    'Position', [0.86986301369863 0.0359712230215827 0.107305936073059 0.12], ...
                    'Callback', {@ApplyNewChanSets, handles}, ...
                    'FontSize', 9, ...
                    'BackgroundColor', [0.831, 0.816, 0.784]);
                
tfh(2).h = uicontrol('Parent', tunFig, 'Units', 'normalized', ...
                    'Style', 'text', ...
                    'Tag', 'ChGrpName', ...
                    'Value', [], 'UserData', slctCh, ...
                    'String', ['channel group: ', num2str(chanSettingGrp(slctCh).grpName)], ...
                    'Position', [0 0.902040816326531 1 0.0938775510204084], ...
                    'Callback', {}, ...
                    'FontSize', 10, 'HorizontalAlignment', 'left', ...
                    'BackgroundColor', [0.9, 0.9, 0.9]);
                
tfh(3).h = uicontrol('Parent', tunFig, 'Units', 'normalized', ...
                    'Style', 'text', ...
                    'Tag', 'LowCutTxt', ...
                    'Value', [], ...
                    'String', 'Hz', ...
                    'Position', [0.148619957537155 0.714285714285714 0.0467091295116773 0.0693877551020408], ...
                    'Callback', {}, ...
                    'FontSize', 9, 'HorizontalAlignment', 'left',...
                    'BackgroundColor', [0.831, 0.816, 0.784]);

tfh(4).h = uicontrol('Parent', tunFig, 'Units','normalized', ...
                    'Style', 'checkbox', ...
                    'Tag', 'LowCutFilter', ...
                    'Value', chanSettingGrp(slctCh).lowCut, ...
                    'String', 'LowCut filter', ...
                    'Position', [0.0114155251141553 0.776567317574517 0.212328767123288 0.1], ...
                    'Callback', {}, ...
                    'FontSize', 9, ...
                    'BackgroundColor', [0.831, 0.816, 0.784]);
                
tfh(5).h = uicontrol('Parent', tunFig, 'Units', 'normalized', ...
                    'Style', 'edit', ...
                    'Tag', 'LowCutFrq', ...
                    'Value', [], ...
                    'String', num2str(chanSettingGrp(slctCh).lowCutFrq), ...
                    'Position', [0.0114155251141553 0.692526794890624 0.136986301369863 0.1], ...
                    'Callback', {}, ...
                    'FontSize', 9, ...
                    'BackgroundColor', [1, 1, 1]);
                
tfh(6).h = uicontrol('Parent', tunFig, 'Units','normalized', ...
                    'Style', 'text', ...
                    'Tag', 'HighCutTxt', ... 
                    'Value', [], ...
                    'String', 'Hz', ...
                    'Position', [0.411889596602972 0.714285714285714 0.0467091295116773 0.0693877551020408], ...
                    'Callback', {}, ...
                    'FontSize', 9, 'HorizontalAlignment', 'left',...
                    'BackgroundColor', [0.831, 0.816, 0.784]);

tfh(7).h = uicontrol('Parent', tunFig, 'Units', 'normalized', ...
                    'Style', 'checkbox', ...
                    'Tag', 'HighCutFilter', ...
                    'Value', chanSettingGrp(slctCh).highCut, ...
                    'String', 'HighCut filter', ...
                    'Position', [0.271689497716895 0.776567317574517 0.212328767123288 0.1], ...
                    'Callback', {}, ...
                    'FontSize', 9, ...
                    'BackgroundColor', [0.831, 0.816, 0.784]);
                
tfh(8).h = uicontrol('Parent', tunFig, 'Units','normalized', ...
                    'Style', 'edit', ...
                    'Tag', 'HighCutFrq', ...
                    'Value', [], ...
                    'String', num2str(chanSettingGrp(slctCh).highCutFrq), ...
                    'Position', [0.271689497716894 0.692526794890624 0.136986301369863 0.1], ...
                    'Callback', {}, ...
                    'FontSize', 9, ...
                    'BackgroundColor', [1, 1, 1]);

tfh(9).h = uicontrol('Parent', tunFig, 'Units', 'normalized', ...
                    'Style', 'checkbox', ...
                    'Tag', 'RecovDC', ...
                    'Value', chanSettingGrp(slctCh).recovDC, ...
                    'String', 'recovDC', ...
                    'Position', [0.0445859872611465 0.567346938775511 0.16135881104034 0.102040816326531], ...
                    'Callback', {}, ...
                    'FontSize', 9, ...
                    'BackgroundColor', [0.831, 0.816, 0.784]);
                
tfh(10).h = uicontrol('Parent', tunFig, 'Units', 'normalized', ...
                    'Style', 'text', ...
                    'Tag', 'LFPscalTxt1', ...
                    'Value', [], ...    
                    'String', 'LFP scale', ...
                    'Position', [0.634819532908705 0.23265306122449 0.176220806794055 0.102040816326531], ...
                    'Callback', {}, ...
                    'FontSize', 10, ...
                    'BackgroundColor', [0.831, 0.816, 0.784]);
                
tfh(11).h = uicontrol('Parent', tunFig, 'Units', 'normalized', ...
                    'Style', 'text', ...
                    'Tag', 'LFPscalTxt2', ...
                    'Value', [], ...
                    'String', 'mV', ...
                    'Position', [0.789808917197452 0.175510204081633 0.0467091295116773 0.0693877551020408], ...
                    'Callback', {}, ...
                    'FontSize', 9, 'HorizontalAlignment', 'left',...
                    'BackgroundColor', [0.831, 0.816, 0.784]);

tfh(12).h = uicontrol('Parent', tunFig, 'Units', 'normalized', ...
                    'Style', 'edit', ...
                    'Tag', 'LFPscale', ...
                    'Value', [], ...
                    'String', num2str(chanSettingGrp(slctCh).lfpScale), ...
                    'Position', [0.658174097664543 0.159183673469388 0.129511677282378 0.102040816326531], ...
                    'Callback', {}, ...
                    'FontSize', 9, ...
                    'BackgroundColor', [1, 1, 1]);
                
tfh(13).h = uicontrol('Parent', tunFig, 'Units', 'normalized', ...
                    'Style', 'pushbutton', ...
                    'Tag', 'ReduceLFP', ...
                    'Value', [], ...
                    'String', '<', ...
                    'Position', [0.651804670912951 0.0489795918367347 0.0658174097664544 0.102040816326531],...
                    'Callback', {@DownsizeAmplScale, tfh(12).h}, ...
                    'BackgroundColor', [0.831, 0.816, 0.784]);

tfh(14).h = uicontrol('Parent', tunFig, 'Units', 'normalized', ...
                    'Style', 'pushbutton', ...
                    'Tag', 'EnlargeLFP', ...
                    'Value', [], ...
                    'String', '>', ...
                    'Position', [0.732484076433121 0.0489795918367347 0.0658174097664544 0.102040816326531], ...
                    'Callback', {@UpsizeAmplScale, tfh(12).h}, ...
                    'BackgroundColor', [0.831, 0.816, 0.784]);

tfh(15).h = uicontrol('Parent', tunFig, 'Units', 'normalized', ...
                    'Style', 'checkbox', ...
                    'Tag', 'CommRef', ...
                    'Value', chanSettingGrp(slctCh).comRef, ...
                    'String', 'Comm.ref.', ...
                    'Position', [0.787685774946921 0.66530612244898 0.182590233545648 0.102040816326531], ...
                    'Callback', [], ...
                    'FontSize', 9, ...
                    'BackgroundColor', [0.831, 0.816, 0.784]);

tfh(16).h = uicontrol('Parent', tunFig, 'Units', 'normalized', ...
                    'Style', 'checkbox',...
                    'Tag', 'PlotCSDmap', ...
                    'Value', chanSettingGrp(slctCh).pltCSDmap, ...
                    'String', 'CSD map', ...
                    'Position', [0.787685774946921 0.551510204081633 0.178343949044586 0.102040816326531], ...
                    'Callback', {@CSDmaporMUAmap, tfh(19).h}, ...
                    'FontSize', 9, ...
                    'BackgroundColor', [0.831, 0.816, 0.784]);

tfh(17).h = uicontrol('Parent', tunFig, 'Units', 'normalized', ...
                    'Style', 'checkbox', ...
                    'Tag', 'DCshift', ...
                    'Value', chanSettingGrp(slctCh).dcShift, ...
                    'String', 'DC shift', ...
                    'Position', [0.787685774946921 0.755102040816327 0.167728237791932 0.102040816326531], ...
                    'Callback', [], ...
                    'FontSize', 9, ...
                    'BackgroundColor', [0.831, 0.816, 0.784]);

tfh(18).h = uicontrol('Parent', tunFig, 'Units', 'normalized', ...
                    'Style', 'checkbox',...
                    'Tag', 'RawSignal', ...
                    'Value', chanSettingGrp(slctCh).rawSignal, ...
                    'String', 'raw signal', ...
                    'Position', [0.222929936305732 0.559183673469388 0.199575371549894 0.102040816326531], ...
                    'Callback', [], ...
                    'FontSize', 9, ...
                    'BackgroundColor', [0.831, 0.816, 0.784]);

tfh(19).h = uicontrol('Parent', tunFig, 'Units', 'normalized', ...
                    'Style', 'checkbox',...
                    'Tag', 'PlotSpkMap', ...
                    'Value', chanSettingGrp(slctCh).pltSpkMap, ...
                    'String', 'Spikes map', ...
                    'Position', [0.787685774946921 0.469387755102041 0.201698513800425 0.102040816326531], ...
                    'Callback', {@CSDmaporMUAmap, tfh(16).h}, ...
                    'FontSize', 9, ...
                    'BackgroundColor', [0.831, 0.816, 0.784]);
set(tfh(16).h, 'Callback', {@CSDmaporMUAmap, tfh(19).h})

tfh(20).h = uicontrol('Parent', tunFig, 'Units', 'normalized', ...
                    'Style', 'checkbox', ...
                    'Tag', 'PlotMUA', ...
                    'Value', chanSettingGrp(slctCh).pltMUA, ...
                    'String', 'MUA', ...
                    'Position', [0.556263269639066 0.755102040816327 0.125265392781316 0.102040816326531], ...
                    'Callback', {}, ...
                    'FontSize', 9, ...
                    'BackgroundColor', [0.831, 0.816, 0.784]);

tfh(21).h = uicontrol('Parent', tunFig, 'Units', 'normalized', ...
                    'Style','checkbox', ...
                    'Tag', 'PlotLFP', ...
                    'Value', chanSettingGrp(slctCh).pltLFP, ...
                    'String', 'LFP', ...
                    'Position', [0.556263269639066 0.575510204081633 0.125265392781316 0.102040816326531], ...
                    'Callback', {}, ... {@MUAorMUAfreq, tfh(22).h}
                    'FontSize', 9, ...
                    'BackgroundColor', [0.831, 0.816, 0.784]);

tfh(22).h = uicontrol('Parent', tunFig, 'Units', 'normalized', ...
                    'Style', 'checkbox', ...
                    'Tag', 'PlotMUAfreq', ...
                    'Value', chanSettingGrp(slctCh).pltMUAfreq, ...
                    'String', 'MUA freq', ...
                    'Position', [0.556263269639066 0.66530612244898 0.182590233545648 0.102040816326531], ...
                    'Callback', {}, ... {@MUAorMUAfreq, tfh(21).h}
                    'FontSize', 9, ...
                    'BackgroundColor', [0.831, 0.816, 0.784]);
%set(tfh(21), 'Callback', {@MUAorMUAfreq, tfh(22).h})

tfh(23).h = uicontrol('Parent', tunFig, 'Units', 'normalized', ...
                    'Style', 'text', ...
                    'Tag', 'MUAbinTxt1', ...
                    'Value', [], ...
                    'String', 'MUA bin', ...
                    'Position', [0.0233545647558386 0.375510204081634 0.114649681528662 0.0700000000000001], ...
                    'Callback', {}, ...
                    'FontSize', 9, ...
                    'BackgroundColor', [0.831, 0.816, 0.784]);

tfh(24).h = uicontrol('Parent', tunFig, 'Units', 'normalized', ...
                    'Style', 'text', ...
                    'Tag', 'MUAbinTxt2', ...
                    'Value', [], ...
                    'String', '(ms)',...
                    'Position', [0.114649681528662 0.289795918367347 0.0658174097664544 0.0775510204081634], ...
                    'Callback', {}, ...
                    'FontSize',9, 'HorizontalAlignment','left', ...
                    'BackgroundColor', [0.831, 0.816, 0.784]);
                
tfh(25).h = uicontrol('Parent', tunFig, 'Units', 'normalized', ...
                    'Style', 'edit', ...
                    'Tag', 'SpkBin', ...
                    'Value', [], ...
                    'String', num2str(chanSettingGrp(slctCh).muaBin), ...
                    'Position', [0.0212314225053079 0.281632653061225 0.0934182590233546 0.102040816326531], ...
                    'Callback', {}, ...
                    'FontSize', 9, ...
                    'BackgroundColor', [1, 1, 1]);

tfh(26).h = uicontrol('Parent', tunFig, 'Units', 'normalized', ...
                    'Style', 'text', ...
                    'Tag', 'PrecisFctTxt', ...
                    'Value', [], ...
                    'String', {  'precision'; 'factor' }, ...
                    'Position', [0.186836518046709 0.379591836734695 0.114649681528662 0.122448979591837], ...
                    'Callback', {}, ...
                    'FontSize', 9, ...
                    'BackgroundColor', [0.831, 0.816, 0.784]);
                
tfh(27).h = uicontrol('Parent', tunFig, 'Units', 'normalized', ...
                    'Style', 'edit', ...
                    'Tag', 'PrecisFct', ...
                    'Value', [], ...
                    'String', num2str(chanSettingGrp(slctCh).precisF),...
                    'Position', [0.199575371549894 0.281632653061225 0.0934182590233546 0.1], ...
                    'Callback', {}, ...
                    'FontSize', 9, ...
                    'BackgroundColor', [1, 1, 1]);

tfh(28).h = uicontrol('Parent', tunFig, 'Units', 'normalized', ...
                    'Style', 'text', ...
                    'Tag', 'MUAscalTxt1', ...
                    'Value', [], ...
                    'String', 'MUA scale', ...
                    'Position', [0.337579617834395 0.310204081632653 0.176220806794055 0.073469387755102], ...
                    'Callback', {}, ...
                    'FontSize', 10, ...
                    'BackgroundColor', [0.831, 0.816, 0.784]);

tfh(29).h = uicontrol('Parent', tunFig, 'Units', 'normalized', ...
                    'Style', 'text', ...
                    'Tag', 'MUAscalTxt2', ...
                    'Value', [], ...
                    'String', 'spik/ms', ...
                    'Position', [0.469214437367304 0.220408163265306 0.106157112526539 0.073469387755102], ...
                    'Callback', {}, ...
                    'FontSize', 9, 'HorizontalAlignment', 'left', ...
                    'BackgroundColor', [0.831, 0.816, 0.784]);
                
tfh(30).h = uicontrol('Parent', tunFig, 'Units', 'normalized', ...
                    'Style','edit', ...
                    'Tag','MUAscale', ...
                    'Value', [], ...
                    'String', num2str(chanSettingGrp(slctCh).muaScale), ...
                    'Position', [0.360934182590234 0.208163265306123 0.106157112526539 0.102040816326531], ...
                    'Callback', {}, ...
                    'FontSize', 9, ...
                    'BackgroundColor', [1, 1, 1]);

tfh(31).h = uicontrol('Parent', tunFig, 'Units', 'normalized', ...
                    'Style', 'pushbutton', ...
                    'Tag','EnlargeMUAfrq',...
                    'Value', [], ...
                    'String', '>', ...
                    'Position', [0.426751592356688 0.0979591836734694 0.0658174097664544 0.102040816326531], ...
                    'Callback', {@UpsizeAmplScale, tfh(30).h}, ...
                    'FontSize', 9, ...
                    'BackgroundColor', [0.831, 0.816, 0.784]);

tfh(32).h = uicontrol('Parent', tunFig, 'Units', 'normalized', ...
                    'Style', 'pushbutton', ...
                    'Tag', 'ReduceMUAfrq',...
                    'Value', [], ...
                    'String', '<',...
                    'Position', [0.346072186836518 0.0979591836734694 0.0658174097664544 0.102040816326531],...
                    'Callback', {@DownsizeAmplScale, tfh(30).h}, ...
                    'FontSize', 9, ...
                    'BackgroundColor', [0.831, 0.816, 0.784]);

tfh(33).h = uicontrol('Parent', tunFig, 'Units', 'normalized', ...
                    'Style', 'text', ...
                    'Tag', 'MUAtrsldTxt1', ...
                    'Value', [], ...
                    'String',{  'MUA'; 'Threshold' }, ...
                    'Position', [0.0828025477707006 0.110204081632653 0.157112526539278 0.126530612244898], ...
                    'Callback', {}, ...
                    'FontSize', 9, ...
                    'BackgroundColor', [0.831, 0.816, 0.784]);

tfh(34).h = uicontrol('Parent', tunFig, 'Units', 'normalized', ...
                    'Style', 'text', ...
                    'Tag', 'MUAtrsldTxt2', ...
                    'Value', [], ...
                    'String', 'std', ...
                    'Position', [0.21656050955414 0.036734693877551 0.059447983014862 0.0653061224489796], ...
                    'Callback', {}, ...
                    'FontSize', 9, 'HorizontalAlignment', 'left', ...
                    'BackgroundColor', [0.831, 0.816, 0.784]);
                
tfh(35).h = uicontrol('Parent', tunFig, 'Units', 'normalized', ...
                    'Style', 'edit', ...
                    'Tag', 'MUAthrshld', ...
                    'Value', [], ...
                    'String', num2str(chanSettingGrp(slctCh).muaTrsld), ...
                    'Position', [0.118895966029724 0.0122448979591837 0.0934182590233546 0.1],...
                    'Callback', {}, ...
                    'FontSize', 9, ...
                    'BackgroundColor', [1, 1, 1]);
                
tfh(36).h = uicontrol('Parent', tunFig, 'Units', 'normalized', ...
                    'Style', 'checkbox', ...
                    'Tag', 'InvertLFP', ...
                    'Value', chanSettingGrp(slctCh).invertLFP, ...
                    'String', 'invert LFP', ...
                    'Position', [0.558386411889597 0.469387755102041 0.169851380042463 0.102040816326531], ...
                    'Callback', {}, ...
                    'FontSize', 9, ...
                    'BackgroundColor', [0.831, 0.816, 0.784]);
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
slctCh = get(tfh(2).h, 'UserData');%original number of selected channel
zavp = getappdata(evi, 'zavp');%additional parameters
lfpOrig = getappdata(evi, 'lfpOrig');%original LFP
lfp = getappdata(evi, 'lfp');%filtered LFP
spksOrig = getappdata(evi, 'spksOrig');
spks = getappdata(evi, 'spks');
lfpVar = getappdata(evi, 'lfpVar');
hd = getappdata(evi, 'hd');
muaprgChang = false;%mua threshold was changed for the group of channels

for ch = 1:length(chanSettingGrp) %run over channels (looking for the group)
    if (chanSettingGrp(ch).grpName == chanSettingGrp(slctCh).grpName) %same group
        nLwCtFr = str2double(get(tfh(5).h, 'String'));%new lowcut frequency
        nHiCtFr = str2double(get(tfh(8).h, 'String'));%new highcut frequency
        newLwCt = ((chanSettingGrp(ch).lowCut ~= get(tfh(4).h, 'Value')) || ... %check or uncheck
                   ((chanSettingGrp(ch).lowCutFrq ~= nLwCtFr) && chanSettingGrp(ch).lowCut && get(tfh(4).h, 'Value')) ...
                  );%lowcut filtering changed
        chanSettingGrp(ch).lowCut = get(tfh(4).h, 'Value');%if do lowcut filtering
        chanSettingGrp(ch).lowCutFrq = nLwCtFr;%lowcut frequency
        
        newHgCt = ((chanSettingGrp(ch).highCut ~= get(tfh(7).h, 'Value')) || ... %check or uncheck
                   ((chanSettingGrp(ch).highCutFrq ~= nHiCtFr) && chanSettingGrp(ch).highCut && get(tfh(7).h, 'Value')) ...
                  );%highcut filtering changed
        chanSettingGrp(ch).highCut = get(tfh(7).h, 'Value');%if do highcut filtering
        chanSettingGrp(ch).highCutFrq = nHiCtFr;%highcut frequency
        
        rercvDC = (chanSettingGrp(ch).recovDC ~= get(tfh(9).h, 'Value'));%change recovery from pseudoDC
        chanSettingGrp(ch).recovDC = get(tfh(9).h, 'Value');%if do DC recovery from pseudoDC
        chanSettingGrp(ch).rawSignal = get(tfh(18).h, 'Value');%if plot raw signal
        
        chanSettingGrp(ch).comRef = get(tfh(15).h, 'Value');%if do common reference subtraction
        chanSettingGrp(ch).dcShift = get(tfh(17).h, 'Value');%if compensate DC
        chanSettingGrp(ch).pltCSDmap = get(tfh(16).h, 'Value');%if plot CSD map
        chanSettingGrp(ch).pltSpkMap = get(tfh(19).h, 'Value');%if plot spikes map
        chanSettingGrp(ch).pltLFP = get(tfh(21).h, 'Value');%if plot LFP traces
        chanSettingGrp(ch).pltMUA = get(tfh(20).h, 'Value');%if plot MUA bars
        chanSettingGrp(ch).pltMUAfreq = get(tfh(22).h, 'Value');%if plot MUA-frequency traces
        chanSettingGrp(ch).lfpScale = str2double(get(tfh(12).h, 'String'));%verticale scale (mV, amplitude)
        newPolar = (chanSettingGrp(ch).invertLFP ~= get(tfh(36).h, 'Value')); %if invert LFP
        chanSettingGrp(ch).invertLFP = get(tfh(36).h, 'Value');%if invert LFP
        
        %= LFP filtration =%
        if (newLwCt || newHgCt || rercvDC) %new filtration parameters
            lfp(:, ch, :) = lfpOrig(:, ch, :);%start from sourse LFP
            if (chanSettingGrp(ch).recovDC) %recovery DC from pseudoDC
                lfp(:, ch, :) = RecovDCfromPseudo(lfp(:, ch, :), zavp.dwnSmplFrq);%recovery
            end
            
            if (chanSettingGrp(ch).lowCut && (chanSettingGrp(ch).lowCutFrq < (zavp.dwnSmplFrq / 2))) %lowcut filtration
                if (chanSettingGrp(ch).lowCutFrq < 2) %small cutting frequency
                    lfp(:, ch, :) = ZavRCfilt(lfp(:, ch, :), chanSettingGrp(ch).lowCutFrq, zavp.dwnSmplFrq, 'high');%RC-filter, lowcut
                else %quite large cutting frequency
                    lfp(:, ch, :) = ZavFilter(lfp(:, ch, :), zavp.dwnSmplFrq, 'high', chanSettingGrp(ch).lowCutFrq, 2);%lowcut
                end
            end
            if (chanSettingGrp(ch).highCut && (chanSettingGrp(ch).highCutFrq < (zavp.dwnSmplFrq / 2))) %highcut filtration
                if (chanSettingGrp(ch).highCutFrq < 2) %small cutting frequency
                    lfp(:, ch, :) = ZavRCfilt(lfp(:, ch, :), chanSettingGrp(ch).highCutFrq, zavp.dwnSmplFrq, 'low');%RC-filter, lowcut
                else
                    lfp(:, ch, :) = ZavFilter(lfp(:, ch, :), zavp.dwnSmplFrq, 'low', chanSettingGrp(ch).highCutFrq, 2);%highcut
                end
            end
        end
        %= end of (LFP filtration) =%
        
        %polarity
        if (newPolar) %if invert LFP
            lfp(:, ch, :) = -1 * lfp(:, ch, :);%invert LFP
        end
        
        chanSettingGrp(ch).muaBin = str2double(get(tfh(25).h, 'String'));%bin step for MUA frequency (ms)
        chanSettingGrp(ch).precisF = str2double(get(tfh(27).h, 'String'));%precision factor (for MUA frequency calculation)
        chanSettingGrp(ch).muaScale = str2double(get(tfh(30).h, 'String'));%MUA-scale (spikes/ms)
    
        %= MUA threshold change =%
        newMUAtrsld = str2double(get(tfh(35).h, 'String'));%new MUA threshold
        if (chanSettingGrp(ch).muaTrsld ~= newMUAtrsld) %MUA threshold was changed
            chanSettingGrp(ch).muaTrsld = newMUAtrsld;%new threshold
            if isfield(spksOrig, 'ampl')
                spks(ch) = spksOrig(ch);%original spikes
                for sw = 1:hd.lActualEpisodes
                    ii = double(spksOrig(ch, sw).ampl) <= (-lfpVar(ch) * newMUAtrsld);
                    spks(ch, sw).tStamp = spksOrig(ch, sw).tStamp(ii);
                    spks(ch, sw).ampl = spksOrig(ch, sw).ampl(ii);
                    if isempty(spks(ch, sw).tStamp)
                        spks(ch, sw).tStamp = zeros(0, 1);
                        spks(ch, sw).ampl = zeros(0, 1);
                    end
                end
                muaprgChang = true;%changes was done
            end
        end
        %= end of (MUA threshold change) =%
    end
end
delete(tunFig)
setappdata(evi, 'chanSettingGrp', chanSettingGrp);%save channels group setting
setappdata(evi, 'lfp', lfp) %refiltered LFP

if muaprgChang %changes was done
    setappdata(evi, 'spks', spks) %accept changes
end

PlotAll(handles);%plot with new parameters



function CSDmaporMUAmap(hObject, eventdata, nObj)
% hObject    handle to Merge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% nObj - checkbox to be uncheked
%
%choose one of map: CSD or MUA
if get(hObject, 'Value') %one of map was choosed
    set(nObj, 'Value', 0)%switch of CSDmap/MUAmap
end



function MUAorMUAfreq(hObject, eventdata, nObj)
% hObject    handle to Merge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% nObj - checkbox to be uncheked
%
%choose one of line type: MUA (verticale lines) or MUA frequency (timecourse)
if get(hObject, 'Value') %one of map was choosed
    set(nObj, 'Value', 0)%switch of CSDmap/MUAmap
end



% --- change LFP scale (downsize) --- %
function DownsizeAmplScale(hObject, eventdata, tObj)
% hObject    handle to Merge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% tObj - target edit-object (change string in)
%
vScl = str2double(get(tObj, 'String')) / 2;%new scale (downsized)
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
set(tObj, 'String', num2str(vScl))%set new value



% --- change LFP scale (upsize) --- %
function UpsizeAmplScale(hObject, eventdata, tObj)
% hObject    handle to Merge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
vScl = str2double(get(tObj, 'String')) * 2;%new scale (upsized)
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
set(tObj, 'String', num2str(vScl))


% --------------------------------------------------------------------
function DeleteEvnt_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to DeleteEvnt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
evi = handles.Eview;%handle of the application window
brst = getappdata(evi, 'brst');
manlDet = getappdata(evi, 'manlDet');%manually detected events

if (~isempty(brst) && isequal(get(handles.ShowBrst, 'State'), 'on')) %bursts exist and switched on
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
elseif (~isempty(manlDet) && get(handles.CutManual, 'Value')) %delete current manyally detected event
    sn = str2double(get(handles.CurrSegm, 'String'));%wanted sweeps or segments
    manlDet(sn) = [];%delete current event
    segms = getappdata(evi, 'segms');%segments of recordation
    segms(sn, :) = [];%delete corresponding segment
    setappdata(evi, 'manlDet', manlDet);%manually detected events
    setappdata(evi, 'segms', segms);%segments of recordation
    set(handles.SegmTxt, 'String', num2str(length(manlDet)));
    set(handles.CutManual, 'String', ['mdet ', num2str(length(manlDet))]);
    PlotAll(handles);%plot with new parameters
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
function SDprmCalc_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to SDprmCalc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

evi = handles.Eview;%main form handle

%calculate parameters of CSD (Cortical Spreading Depression)
brst = getappdata(evi, 'brst');
if ~isempty(brst) %bursts exist
    pth = getappdata(evi, 'pth');%directory
    flNm = getappdata(evi, 'flNm');%filename
    
    lfp = getappdata(evi, 'lfp');
    lfp = double(lfp);
    spks = getappdata(evi, 'spks');
    %lfpVar = getappdata(evi, 'lfpVar');
    hd = getappdata(evi, 'hd');
    zavp = getappdata(evi, 'zavp');
    prm.Fs = 1 / zavp.siS;%sampling frequency for raw data
    
    jj = vertcat(brst(:).t);
    [~, ii] = sort(jj(:, 1));
    brst = brst(ii);%sorted bursts
    sdprms(1:length(brst)) = struct('p', NaN(16, 21));%parameters of all SDs
    for ch = unique(horzcat(brst(:).ch)) %run over channels
        fLFP = ZavFilter(lfp(:, ch, 1), prm.Fs, 'low', 5, 2);
        fLFP = fLFP - median(fLFP);
        whtLFP = ZavRCfilt(fLFP, 0.001, prm.Fs, 'high');
        for t = 1:length(brst) %run over bursts
            if ismember(ch, brst(t).ch) %SD on current channel
                jj = brst(t).t(1):brst(t).t(2);
                trac = smooth(whtLFP(jj), 2001);
                trac = smooth(ZavMultiDiff(trac, 200), 1901);%smoothed derivation
                [~, p1] = min(trac);%point of maximal slope of SD front
                try
                    sdprms(t).p(ch, :) = ZavGetSDprms(whtLFP(jj), [jj(1), p1], lfp(jj, ch, 1), spks(ch).tStamp, 0, ch);
                    %                    ZavGetSDprms(trac, sgSt1, orgTrc, spkTS, plotFlg, rCh)
                    sdprms(t).p(ch, [1, 2, 4:6, 10, 12]) = sdprms(t).p(ch, [1, 2, 4:6, 10, 12]) * 1e-3;%convert to seconds
                catch
                    sdprms(t).p(ch, :) = -Inf;
                end
            end
        end
    end
    
    %export parameter of SDs to xls-file
    % 1 - time of SD onset; 2 - time of SD negative peak; 3 - amplitude of SD (uV)
    % 4 - time of SD end; 5 - width of SD on half height (ms); 6 - length of SD (full length)
    % 7 - MUA during SD (units); 8 - MUA frequency during SD (1/s); 9 - channel;
    % 10 - time of SD front slope (ms from segment begin); 11 - SD front slope (uV/ms)
    % 12 - time of afterhyperpolarization (AHP); 13 - amplitude of AHP;
    % 14 - T1/2-surface; 15 - T2/2 ("double" T1/2-surface); 16 - negative 3s-plato level
    
    paramToXLS = cell(length(sdprms) * 1000, 21);
    paramToXLS(1, :) = {'file', 'SD #', 'start time (s)', 'channel', ... 1-4
                        'mua backgrnd (1/s)', 'mua during SD (units)', 'mua freq during SD (1/s)', ... 5-7
                        'SD onset (s)', 'SD neg peak (s)', 'SD end (s)', 'time of SD front maximal slope (s)', 'time of AHP (s)', ... 8-12
                        'SD ampl (uV)', 'width of SD on half height (s)', 'full length of SD (s)', ... 13-15
                        'SD front slope (uV/ms)', 'AHP ampl (uV)', 'T1/2-surface (V*s)', 'normalized slope (1/s)', 'T2/2 (V*s)', 'SD plato, uV'};%16-21

    m = 2;%lines counter
    for z = 1:length(sdprms) %run over SDs
        prm = sdprms(z).p;%current SD
        for nch = 1:length(prm(:, 9)) %run over channels
            paramToXLS{m, 1} = zavp.file;%filename
            paramToXLS{m, 2} = z;%SD local # (local number of SD)
            paramToXLS{m, 3} = min(prm(:, 10)) + hd.recTime(1);%start time (s)
            paramToXLS{m, 4} = prm(nch, 9);%channel

            paramToXLS{m, 5} = numel(spks(nch).tStamp) / diff(hd.recTime);%muaBF(nch);%mua backgrnd (1/s)
            paramToXLS{m, 6} = prm(nch, 7);%mua during SD (units)
            paramToXLS{m, 7} = prm(nch, 8);%mua during SD freq (1/s)

            paramToXLS{m, 8} = prm(nch, 1);%SD onset (s)
            paramToXLS{m, 9} = prm(nch, 2);%SD neg peak (s)
            paramToXLS{m, 10} = prm(nch, 4);%SD end (s)
            paramToXLS{m, 11} = prm(nch, 10);%time of SD front maximal slope (s)
            paramToXLS{m, 12} = prm(nch, 12);%time of AHP (s)

            paramToXLS{m, 13} = prm(nch, 3);%SD ampl (uV)
            paramToXLS{m, 14} = prm(nch, 5);%width of SD on half height (s)
            paramToXLS{m, 15} = prm(nch, 6);%full length of SD (s)

            paramToXLS{m, 16} = prm(nch, 11);%SD front slope (uV/ms)
            paramToXLS{m, 17} = prm(nch, 13);%AHP ampl (uV)
            paramToXLS{m, 18} = 1e-6 * prm(nch, 14);%T1/2-surface (V*s)
            paramToXLS{m, 19} = 1e3 * prm(nch, 11) / prm(nch, 3);%normalized lope (1/s)
            paramToXLS{m, 20} = 1e-6 * prm(nch, 15);%T2/2-surface (V*s)
            paramToXLS{m, 21} = prm(nch, 16);%SD negative plato (uV, 3s-plato level)

            m = m + 1;%next line
        end
    end
    paramToXLS(m:end, :) = [];%delete empty
    xlswrite([pth, flNm, '_SDprms.xlsx'], paramToXLS, 1);
end

% --------------------------------------------------------------------
function ManualDetect_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to ManualDetect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
%manual detection
evi = handles.Eview;%main form handle
lfpLn = getappdata(evi, 'lfpLn');%lfp-line handles
if ~isempty(lfpLn) 
    manDetPrmH = handles.ManDetPrm;%table of manual detection parameters
    if isequal(get(manDetPrmH, 'Visible'), 'on') %accept changes
        set(manDetPrmH, 'Visible', 'off')
        set(handles.AcceptManDet, 'Visible', 'off')
        set(evi, 'WindowButtonDownFcn', {@Eview_WindowButtonDownFcn, handles})%reaction on button down
        set(evi, 'WindowButtonMotionFcn', {@Eview_WindowButtonMotionFcn, handles})%reaction on button motion
        set(evi, 'WindowButtonUpFcn', {@Eview_WindowButtonUpFcn, handles})%reaction on button up
        PlotAll(handles);%plot with new parameters
    else %show table of manual detection parameters
        set(manDetPrmH, 'Visible', 'on', 'Position', [0.23714, 0.80854, 0.33, 0.19008])
        set(handles.AcceptManDet, 'Visible', 'on', 'Position', [0.27898, 0.82058, 0.06015, 0.0303])

        set(evi, 'WindowButtonDownFcn', {})%no reaction on button down
        set(evi, 'WindowButtonMotionFcn', {})%no reaction on button motion
        set(evi, 'WindowButtonUpFcn', {})%no reaction on button up
    end
end



% --------------------------------------------------------------------
function ExportManDet_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to ExportManDet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
%export parameters of manually detected events
evi = handles.Eview;%main form handle
pth = getappdata(evi, 'pth');%directory
flNm = getappdata(evi, 'flNm');%filename
zavp = getappdata(evi, 'zavp');%additional parameters
hd = getappdata(evi, 'hd');%header

manlDet = getappdata(evi, 'manlDet');%manually detected events
if ~isempty(manlDet) %events was detected
    isSwps = isfield(hd, 'sweepStartInPts');%sweep mode
    lfp = getappdata(evi, 'lfp');%LFP
    whtLFP = ZavRCfilt(lfp(:, [manlDet(1).ch, manlDet(1).subCh], :), 1, zavp.dwnSmplFrq, 'high');%filtering
    paramToXLS = cell(length(manlDet) + 1, 15);%parameters matrix
    paramToXLS(1, :) = {'Event #', 'Channel', 'Time (s)', 'Ampl (uV)', 'Left peak (ms)', 'Left ampl (uV)', 'Right peak (ms)', 'Right ampl (uV)', ...
                        'SubChannel', 'dT (s)', 'SumAmpl (uV)', 'SubLeft peak (ms)', 'SubLeft ampl (uV)', 'SubRight peak (ms)', 'SubRight ampl (uV)'};%

    for m = 1:length(manlDet) %run over manually detected events
        sw = manlDet(m).sw;%sweep number
        p1 = ZavFindMins(-whtLFP(:, 1, sw));%local maxima
        p2 = ZavFindMins(-whtLFP(:, 1, sw));%local maxima
        paramToXLS{m + 1, 1} = m;%event number
        paramToXLS{m + 1, 2} = manlDet(m).ch;%channel where event was detected
        p0 = 0;%sweep start (seconds from beginning of recording)
        if isSwps %sweep start exist
            p0 = hd.sweepStartInPts(sw) * hd.si / 1e6;%sweep start (seconds from beginning of recording)
        end
        paramToXLS{m + 1, 3} = (manlDet(m).t / 1e3) + p0;%time of event (s from record begin)
        manlDet(m).t = round(manlDet(m).t);
        paramToXLS{m + 1, 4} = whtLFP(manlDet(m).t, 1, sw);%absolute amplitude (uV)
        
        t = find(p1 < manlDet(m).t, 1, 'last');%left maximum
        if ~isempty(t)
            paramToXLS{m + 1, 5} = p1(t) - manlDet(m).t;%left maximum
            paramToXLS{m + 1, 6} = whtLFP(p1(t), 1, sw);%left amplitude (uV)
        end
        t = find(p1 > manlDet(m).t, 1, 'first');%right maximum
        if ~isempty(t)
            paramToXLS{m + 1, 7} = p1(t) - manlDet(m).t;%right maximum
            paramToXLS{m + 1, 8} = whtLFP(p1(t), 1, sw);%right amplitude (uV)
        end
        
        %= subchannel data =%
        if (~isempty(manlDet(m).subCh) && ~isempty(manlDet(m).subT)) %subchannel was treated
            paramToXLS{m + 1, 9} = manlDet(m).subCh;%subchannel
            paramToXLS{m + 1, 10} = round(manlDet(m).subT - manlDet(m).t);%time delta of event on subchannel (ms from record begin)
            manlDet(m).subT = round(manlDet(m).subT);
            paramToXLS{m + 1, 11} = whtLFP(manlDet(m).subT, 2, sw);%absolute amplitude on subchannel (uV)

            t = find(p2 < manlDet(m).subT, 1, 'last');%left maximum on subchannel
            if ~isempty(t)
                paramToXLS{m + 1, 12} = p2(t) - manlDet(m).subT;%left maximum on subchannel
                paramToXLS{m + 1, 13} = whtLFP(p2(t), 2, sw);%left amplitude on subchannel (uV)
            end
            t = find(p2 > manlDet(m).subT, 1, 'first');%right maximum on subchannel
            if ~isempty(t)
                paramToXLS{m + 1, 14} = p2(t) - manlDet(m).subT;%right maximum on subchannel
                paramToXLS{m + 1, 15} = whtLFP(p2(t), 2, sw);%right amplitude on subchannel (uV)
            end
        end
    end
    paramToXLS((length(manlDet) + 1):end, :) = [];%delete empty
    xlswrite([pth, flNm, '_ManlDetEvntPrm.xlsx'], paramToXLS, 1);
end



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
[imgFlNm, pth] = uiputfile({'*.emf', 'EMF'; '*.eps', 'EPS'; '*.jpg', 'JPG'}, 'Save as', [pth, flNm]);%save image dialog
z = find(imgFlNm == '.', 1, 'last');
imtyp = imgFlNm((z + 1):end);%wanted image type

if isequal(imtyp, 'emf') %EMF
    saveas(evi, [pth, imgFlNm], imtyp)%save image
elseif isequal(imtyp, 'jpg') %JPG
    imgM = getframe(evi);%get objects of picture
    imwrite(imgM.cdata, [pth, imgFlNm], imtyp)%save image
elseif isequal(imtyp, 'eps') %EMF
    saveas(evi, [pth, imgFlNm], imtyp)%save image
end



% --- Executes on button press in SelectChann.
function SelectChann_Callback(hObject, eventdata, handles)
% hObject    handle to SelectChann (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
%channels visibility settings
evi = handles.Eview;%main form handle
lfpLn = getappdata(evi, 'lfpLn');%lfp-line handles
if ~isempty(lfpLn) 
    chanSelecH = handles.ChannSelector;%table of channels selection parameters
    if isequal(get(chanSelecH, 'Visible'), 'on') %accept changes
        set(chanSelecH, 'Visible', 'off')
        set(handles.SelectAllCh, 'Visible', 'off')
        set(handles.InvertSelection, 'Visible', 'off')
        set(evi, 'WindowButtonDownFcn', {@Eview_WindowButtonDownFcn, handles})%reaction on button down
        set(evi, 'WindowButtonMotionFcn', {@Eview_WindowButtonMotionFcn, handles})%reaction on button motion
        set(evi, 'WindowButtonUpFcn', {@Eview_WindowButtonUpFcn, handles})%reaction on button up
        PlotAll(handles);%plot with new parameters
    else %show channels selector
        set(chanSelecH, 'Visible', 'on')
        set(handles.SelectAllCh, 'Visible', 'on')
        set(handles.InvertSelection, 'Visible', 'on')

        set(evi, 'WindowButtonDownFcn', {})%no reaction on button down
        set(evi, 'WindowButtonMotionFcn', {})%no reaction on button motion
        set(evi, 'WindowButtonUpFcn', {})%no reaction on button up
    end
end



% --- Executes on button press in SelectAllCh.
function SelectAllCh_Callback(hObject, eventdata, handles)
% hObject    handle to SelectAllCh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

chanSelecH = handles.ChannSelector;%table of channels selections
chTabData = get(chanSelecH, 'Data');%channel selection table
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
chTabData = get(chanSelecH, 'Data');%channel selection table
for t = 1:size(chTabData, 1) %run over rows of table (channels)
    chTabData{t, 3} = ~chTabData{t, 3};%invert the selection
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
evi = handles.Eview;%main form handle
lfpLn = getappdata(evi, 'lfpLn');%lfp-line handles
spkDnLn = getappdata(evi, 'spkDnLn');%spike density line handles

tm_p = [];
lfp = [];
tm_s = [];
spkD = [];

if ~isempty(lfpLn) %LFP visible
    tm_p = get(lfpLn(1), 'XData');%lfp time vector
    lfp = NaN(length(tm_p), length(lfpLn));
    for t = 1:length(lfpLn) %run over LFP-traces
        if isequal(get(lfpLn(t), 'Visible'), 'on')
            dCoeffP = get(lfpLn(t), 'UserData');
            ch = dCoeffP{5};%actual channel
            lfp(:, ch) = (get(lfpLn(t), 'YData') - dCoeffP{2} - dCoeffP{4}) * dCoeffP{1}(3);
        end
    end
    tm_p = tm_p * dCoeffP{3};%lfp time vector in milliseconds
end

if ~isempty(spkDnLn) %spikes visible
    tm_s = get(spkDnLn(1), 'Xdata');%spikes time vector
    spkD = NaN(length(tm_s), length(spkDnLn));%spikes density
    for t = 1:length(spkDnLn) %run over MUAdense-traces
        if isequal(get(spkDnLn(t), 'Visible'), 'on')
            dCoeffP = get(spkDnLn(t), 'UserData');
            ch = dCoeffP{5};%actual channel
            spkD(:, ch) = (get(spkDnLn(t), 'YData') - dCoeffP{2} - dCoeffP{4}) * dCoeffP{1}(3);
        end
    end
    tm_s = tm_s * dCoeffP{3};%spikes time vector in milliseconds
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
eval([spkName, ' = [tm_s'', spkD];'])%spikes values

%%= save to warkspace =%%
eval(['assignin(''base'', ''', lfpName, ''', ', lfpName, ')'])
eval(['assignin(''base'', ''', spkName, ''', ', spkName, ')'])
%%= end of (save to warkspace) =%%

% %%= save to mat-file =%%
% if exist('eview_screen.mat', 'file') %storage exist
%     save('eview_screen.mat', lfpName, '-append')
% else %no storage
%     save('eview_screen.mat', lfpName)
% end
% if eval(['~isempty(', spkName, ')']) %spikes are shown
%     save('eview_screen.mat', spkName, '-append')
% end
% %%= end of (save to mat-file) =%%




% --- Executes when selected cell(s) is changed in ChannSelector.
function ChannSelector_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to ChannSelector (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)

chanSelecH = handles.ChannSelector;%table of channels selections
ii = unique(eventdata.Indices(:, 1));%indices of selected rows
set(chanSelecH, 'UserData', eventdata.Indices)%save current selection
if ((length(ii) > 1) && all(diff(ii) == 1) && all(eventdata.Indices(:, 2) == 3)) %change channels selection
    chTabData = get(chanSelecH, 'Data');%read current table content
    jj = horzcat(chTabData{ii, 3});%selection falgs
    trg = ~(sum(jj) >= sum(~jj));%target value (prevailing selection)
    for z = ii' %run over selected raws
        chTabData{z, 3} = trg;%change value of selected cells (set prevailing value)
    end
    set(chanSelecH, 'Data', chTabData)%set table data
    
%     jscrollpane = javaObjectEDT(get(chanSelecH, 'UserData'));%java-handle of scrollpane
%     viewport = javaObjectEDT(jscrollpane.getViewport);%jave-handle
%     jtable = javaObjectEDT(viewport.getView);%jave-handle of table content
%     jtable.scrollRowToVisible(ii(end));%make selected rows visible again
end

if ((length(ii) > 1) && all(diff(ii) == 1) && all(eventdata.Indices(:, 2) == 4)) %change channels groups
    evi = handles.Eview;%handle of the application window
    chTabData = get(chanSelecH, 'Data');%read current table content
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
    set(chanSelecH, 'Data', chTabData)%set table data
    setappdata(evi, 'chanSettingGrp', chanSettingGrp);%save channels group setting
end


function SeparateChannel(hObject, eventdata, handles)
%creat individual group for the selected channel
%
evi = handles.Eview;%handle of the application window
chanSelecH = handles.ChannSelector;%table of channels selections
chTabData = get(chanSelecH, 'Data');%read current table content
chanSettingGrp = getappdata(evi, 'chanSettingGrp');%load channels group setting
ii = get(chanSelecH, 'UserData');%load last selection
ii = unique(ii(:, 1));%indices of selected rows
m = min(setdiff(1:size(chTabData, 1), horzcat(chTabData{:, 4})));%number of new unique group to be add

for z = ii' %run over selected rows (channels)
    chanSettingGrp(z).grpName = m;%name of group for all selected channels
    chTabData{z, 4} = chanSettingGrp(z).grpName;%copy group name
end
set(chanSelecH, 'Data', chTabData)%set table data
setappdata(evi, 'chanSettingGrp', chanSettingGrp);%save channels group setting
    


function CSDcMapMax_Callback(hObject, eventdata, handles)
% hObject    handle to CSDcMapMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%change color axis
set(handles.CSDcMapMin, 'String', ['-', get(handles.CSDcMapMax, 'String')]) %set minimum automatically
caxis(handles.LFP_ax, [str2double(get(handles.CSDcMapMin, 'String')), str2double(get(handles.CSDcMapMax, 'String'))])%colorbar range



function SpkCMap_Callback(hObject, eventdata, handles)
% hObject    handle to SpkCMapMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%change color axis
caxis(handles.LFP_ax, [str2double(get(handles.SpkCMapMin, 'String')), str2double(get(handles.SpkCMapMax, 'String'))])%colorbar range



function data = RemovNaN(data)
%data = RemovNaN(data)
%exclude NaN and Inf from LFP with interpolation
%
%INPUTS
%data - original signal
%
%OUTPUS
%data - signal without NaN and Inf
%
data(isinf(data)) = NaN;
for sw = 1:size(data, 3) %run over sweeps
    for ch = 1:size(data, 2) %run over channels
        ii = (1:size(data, 1)) / 1e4;
        ii(isnan(data(:, ch, sw))) = -10;
        [pntsL, pntsR] = ZavFindFlats(ii, 1e-7);%flat segments
        for t = 1:length(pntsL) %run over NaN-segments
            jj1 = pntsL(t) + (-50:-1);
            jj1((jj1 < 1) | (jj1 > size(data, 1))) = [];
            jj2 = pntsR(t) + (1:50);
            jj2((jj2 < 1) | (jj2 > size(data, 1))) = [];
            data(pntsL(t):pntsR(t), ch, sw) = interp1([pntsL(t); pntsR(t)], [nanmedian(data(jj1, ch, sw)), nanmedian(data(jj2, ch, sw))], (pntsL(t):pntsR(t))');
        end
        ii = find(isnan(data(:, ch, sw)));%solitary inf-nan points
        for t = 1:length(ii) %run over solitary inf-nan points
            jj1 = ii(t) + (-50:-1);
            jj1((jj1 < 1) | (jj1 > size(data, 1))) = [];
            jj2 = ii(t) + (1:50);
            jj2((jj2 < 1) | (jj2 > size(data, 1))) = [];
            data(ii(t), ch, sw) = interp1([ii(t) - 1; ii(t) + 1], [nanmedian(data(jj1, ch, sw)), nanmedian(data(jj2, ch, sw))], ii(t));
        end
    end
end


% --- Executes on button press in ShowMarkers.
function ShowMarkers_Callback(hObject, eventdata, handles)
% hObject    handle to ShowMarkers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

PlotAll(handles);
