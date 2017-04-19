% Main SISTool file
%
% Coded by: Ali Haddad

function varargout = SISTool(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SISTool_OpeningFcn, ...
                   'gui_OutputFcn',  @SISTool_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
   gui_State.gui_Callback = str2func(varargin{1});
end
if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
end

function SISTool_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);
end

function varargout = SISTool_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;
end

function Load_Callback(hObject, eventdata, handles)
[fname,pname] = uigetfile('*.mat');
if ischar(fname) && ischar(pname)
    load([pname,fname],'EEG')
    if exist('EEG','var')
        handles = load_EEG(handles,EEG,[pname,fname]);
    else
        errordlg('EEG data structure does not exist in the selected file!','Loading Error')
    end
end
guidata(hObject, handles);
end

function Save_Callback(hObject, eventdata, handles)
if isfield(handles,'EEG')
    EEG           = handles.EEG;
    [fname,pname] = uiputfile('*.mat');
    if ischar(fname) && ischar(pname)
        save([pname,fname],'EEG')
        handles = set_fname(handles,[pname,fname]);
    end
else
    errordlg('No loaded EEG data structure!','Saving Error')
end
guidata(hObject, handles);
end

function Import_Callback(hObject, eventdata, handles)
ws = evalin('base','who');
if ~isempty(ws)
    arrname = '';
    d       = dialog('Position',[300 300 250 160]);
    uicontrol('Parent',d,'Style','text','Position',[20 80 210 50],'String','Select MATLAB workspace variable to import:');
    uicontrol('Parent',d,'Position',[89 20 70 25],'String','cancel','Callback','delete(gcf)');
    uicontrol('Parent',d,'Style','popup','Position',[75 70 100 25],'String',ws,'Callback',@popup_callback);
    uiwait(d);
    if ~isempty(arrname)
        EEG = evalin('base',arrname);
        if isnumeric(EEG)
            handles = import_EEG(handles,EEG,arrname);
        else
            handles = load_EEG(handles,EEG,arrname);
        end
    end
else
    errordlg('No variables to import from MATLAB workspace!','Import Error')
end
guidata(hObject, handles);

    function popup_callback(popup,event)
        strings = get(popup,'String');
        arrname = strings{get(popup,'Value')};
        delete(gcf)
    end

end

function Export_Callback(hObject, eventdata, handles)
if isfield(handles,'EEG')
    arrname = {'0'};
    while ~isempty(arrname) && (isempty(arrname{1}) || ~isletter(arrname{1}(1)))
        arrname = inputdlg('Export to MATLAB workspace as:');
    end
    if ~isempty(arrname)
        assignin('base',arrname{1},handles.EEG);
    end
else
    errordlg('No loaded EEG data structure!','Export Error')
end
guidata(hObject, handles);
end

function Segment_Callback(hObject, eventdata, handles)
cond1 = ~isempty(get(handles.start_time,'String')) && ~isempty(get(handles.srate,'String'));
cond2 = ~isempty(get(handles.Wr,'String')) && ~isempty(get(handles.Wd,'String'));
if isfield(handles,'EEG') && cond1 && cond2
    EEGdata    = handles.EEG.data;
    start_time = str2double(get(handles.start_time,'String'));
    srate      = str2double(get(handles.srate,'String'));
    Wr         = round(str2double(get(handles.Wr,'String')) * srate / 1000);
    Wd         = round(str2double(get(handles.Wd,'String')) * srate / 1000);
    if ~isempty(get(handles.Ws,'String'))
        Ws = round(str2double(get(handles.Ws,'String')) * srate / 1000);
    else
        Ws = Wr;
    end
    if ~isempty(get(handles.Wp,'String'))
        Wp = round(str2double(get(handles.Wp,'String')) * srate / 1000);
    else
        Wp = 1;
    end
    if ~isempty(get(handles.Wv,'String'))
        Wv = round(str2double(get(handles.Wv,'String')) * srate / 1000);
    else
        Wv = 0;
    end
    if str2double(get(handles.ntrials,'String')) == 1
        segpnts{1} = SISegmentation(EEGdata,[Wr,Ws,Wp,Wv],Wd);
    else
        segpnts = SISegmentation(EEGdata,[Wr,Ws,Wp,Wv],Wd);
    end
    for trial = 1:length(segpnts)
        segpnts{trial} = start_time + (segpnts{trial}-1) * 1000 / srate;
    end
    handles = set_segpnts(handles,segpnts);
else
    errordlg('EEG data, info entries, or required parameters are missing!','Segmentation Error')
end
guidata(hObject, handles);
end

function Plot_Callback(hObject, eventdata, handles)
if isfield(handles,'EEG') && strcmp(get(handles.segindicator,'String'),"true") ...
        && isfield(handles.EEG,'time')
    EEGdata    = handles.EEG.data;
    ntrials    = str2double(get(handles.ntrials,'String'));
    segpnts    = handles.EEG.segpnts;
    time       = handles.EEG.time;
    if ntrials == 1
        trial = 1;
    else
        trial = 0;
        while ~isempty(trial) && ((trial~=round(trial)) || (trial<=0) || (trial>ntrials))
            trial = str2double(inputdlg(['Select a trial [1 - ',num2str(ntrials),']']));
        end
    end
    if ~isempty(trial)
        figure, plot(time,EEGdata(:,:,trial)')
        pnts = [segpnts{trial},time(end)];
        set(gca,'xlim',[time(1),time(end)],'XTick',pnts,'XTickLabel',num2str(pnts'))
        set(gca,'XGrid','on','FontSize',14,'Fontweight','bold','GridAlpha',1)
        xlabel('segment boundaries [msec]','Fontsize',20,'Fontweight','bold')
        strings = get(handles.Amplabel,'String');
        ylabel(['amplitude [',strings{get(handles.Amplabel,'Value')},']'],'Fontsize',20,'Fontweight','bold')
    end
else
    errordlg('EEG data, segmentation structure, time structure, or info entries are missing!','Plot Error')
end
guidata(hObject, handles);
end

function Amplabel_Callback(hObject, eventdata, handles)
end

function Amplabel_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function start_time_Callback(hObject, eventdata, handles)
start_time = str2double(get(handles.start_time,'String'));
if ~isnan(start_time)
    handles = set_start_time(handles,start_time);
else
    handles = set_start_time(handles,[]);
end
guidata(hObject, handles);
end

function start_time_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function srate_Callback(hObject, eventdata, handles)
srate = str2double(get(handles.srate,'String'));
if ~isnan(srate)
    handles = set_srate(handles,srate);
    handles = update_windows(handles);
else
    handles = set_srate(handles,[]);
end
guidata(hObject, handles);
end

function srate_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function Wr_Callback(hObject, eventdata, handles)
Wr = str2double(get(handles.Wr,'String'));
if ~isnan(Wr)
    handles = set_Wr(handles,Wr);
    handles = update_windows(handles);
else
    handles = set_Wr(handles,[]);
end
guidata(hObject, handles);
end

function Wr_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function Wd_Callback(hObject, eventdata, handles)
Wd = str2double(get(handles.Wd,'String'));
if ~isnan(Wd)
    handles = set_Wd(handles,Wd);
    handles = update_windows(handles);
else
    handles = set_Wd(handles,[]);
end
guidata(hObject, handles);
end

function Wd_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function Ws_Callback(hObject, eventdata, handles)
Ws = str2double(get(handles.Ws,'String'));
if ~isnan(Ws)
    handles = set_Ws(handles,Ws);
    handles = update_windows(handles);
else
    handles = set_Ws(handles,[]);
end
guidata(hObject, handles);
end

function Ws_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function Wp_Callback(hObject, eventdata, handles)
Wp = str2double(get(handles.Wp,'String'));
if ~isnan(Wp)
    handles = set_Wp(handles,Wp);
    handles = update_windows(handles);
else
    handles = set_Wp(handles,[]);
end
guidata(hObject, handles);
end

function Wp_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function Wv_Callback(hObject, eventdata, handles)
Wv = str2double(get(handles.Wv,'String'));
if ~isnan(Wv)
    handles = set_Wv(handles,Wv);
    handles = update_windows(handles);
else
    handles = set_Wv(handles,[]);
end
guidata(hObject, handles);
end

function Wv_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function handles = load_EEG(handles,EEG,fname)
if isstruct(EEG) && isfield(EEG,'data')
    [handles,succ] = import_EEG(handles,EEG.data,fname);
    if succ
        if isfield(EEG,'start_time') && isnumeric(EEG.start_time)
            handles = set_start_time(handles,EEG.start_time);
        end
        if isfield(EEG,'srate') && isnumeric(EEG.srate)
            handles = set_srate(handles,EEG.srate);
        end
        if isfield(EEG,'Wr') && isnumeric(EEG.Wr)
            handles = set_Wr(handles,EEG.Wr);
        end
        if isfield(EEG,'Wd') && isnumeric(EEG.Wd)
            handles = set_Wd(handles,EEG.Wd);
        end
        if isfield(EEG,'Ws') && isnumeric(EEG.Ws)
            handles = set_Ws(handles,EEG.Ws);
        end
        if isfield(EEG,'Wp') && isnumeric(EEG.Wp)
            handles = set_Wp(handles,EEG.Wp);
        end
        if isfield(EEG,'Wv') && isnumeric(EEG.Wv)
            handles = set_Wv(handles,EEG.Wv);
        end
        handles = update_windows(handles);
        if isfield(EEG,'segpnts')
            handles = set_segpnts(handles,EEG.segpnts);
        else
            handles = set_segindicator(handles,false);
        end
    end
else
    errordlg('Unrecognized EEG data structure!','Loading Error')
end
end

function [handles,succ] = import_EEG(handles,EEGdata,arrname)
if isnumeric(EEGdata)
    handles                     = reset_parameters(handles);
    [nchannels,nframes,ntrials] = size(EEGdata);
    handles.EEG.data            = EEGdata;
    handles                     = set_fname(handles,arrname);
    handles                     = set_nchannels(handles,nchannels);
    handles                     = set_nframes(handles,nframes);
    handles                     = set_ntrials(handles,ntrials);
    succ                        = true;
else
    succ = false;
    errordlg('Unrecognized EEG data structure!','Loading Error')
end
end

function handles = set_fname(handles,fname)
set(handles.fname,'String',fname)
end

function handles = set_nchannels(handles,nchannels)
set(handles.nchannels,'String',num2str(nchannels))
end

function handles = set_nframes(handles,nframes)
set(handles.nframes,'String',num2str(nframes))
end

function handles = set_ntrials(handles,ntrials)
set(handles.ntrials,'String',num2str(ntrials))
end

function handles = set_start_time(handles,start_time)
set(handles.start_time,'String',num2str(start_time))
end

function handles = set_srate(handles,srate)
if ~isempty(srate)
    if srate > 0
        set(handles.srate,'String',num2str(srate))
    else
        set(handles.srate,'String','')
        errordlg('Sampling rate has to be positive','Parameter Error')
    end
else
    set(handles.srate,'String','')
end
end

function handles = set_Wr(handles,Wr)
if ~isempty(Wr)
    cond = isempty(get(handles.Wv,'String')) || (str2double(get(handles.Wv,'String'))<Wr);
    if (Wr>0) && cond
        set(handles.Wr,'String',num2str(Wr))
    else
        set(handles.Wr,'String','')
        errordlg('Reference window has to be positive and greater than overlap (if defined)','Parameter Error')
    end
else
    set(handles.Wr,'String','')
end
end

function handles = set_Wd(handles,Wd)
if ~isempty(Wd)
    if Wd > 0
        set(handles.Wd,'String',num2str(Wd))
    else
        set(handles.Wd,'String','')
        errordlg('Decision window has to be positive','Parameter Error')
    end
else
    set(handles.Wd,'String','')
end
end

function handles = set_Ws(handles,Ws)
if ~isempty(Ws)
    cond = isempty(get(handles.Wv,'String')) || (str2double(get(handles.Wv,'String'))<Ws);
    if (Ws>0) && cond
        set(handles.Ws,'String',num2str(Ws))
    else
        set(handles.Ws,'String','')
        errordlg('Sliding window has to be positive and greater than overlap (if defined)','Parameter Error')
    end
else
    set(handles.Ws,'String','')
end
end

function handles = set_Wp(handles,Wp)
if ~isempty(Wp)
    if Wp > 0
        set(handles.Wp,'String',num2str(Wp))
    else
        set(handles.Wp,'String','')
        errordlg('Step has to be positive','Parameter Error')
    end
else
    set(handles.Wp,'String','')
end
end

function handles = set_Wv(handles,Wv)
if ~isempty(Wv)
    cond1 = isempty(get(handles.Ws,'String')) || (Wv<str2double(get(handles.Ws,'String')));
    cond2 = isempty(get(handles.Wr,'String')) || (Wv<str2double(get(handles.Wr,'String')));
    if (Wv>=0) && cond1 && cond2
        set(handles.Wv,'String',num2str(Wv))
    else
        set(handles.Wv,'String','')
        errordlg('Overlap has to be non-negative and less than reference and sliding windows (if defined)','Parameter Error')
    end
else
    set(handles.Wv,'String','')
end
end

function handles = update_windows(handles)
if ~isempty(get(handles.srate,'String'))
    srate = str2double(get(handles.srate,'String'));
    if ~isempty(get(handles.Wr,'String'))
        Wr      = round(str2double(get(handles.Wr,'String')) * srate / 1000) * 1000 / srate;
        handles = set_Wr(handles,Wr);
    end
    if ~isempty(get(handles.Wd,'String'))
        Wd      = round(str2double(get(handles.Wd,'String')) * srate / 1000) * 1000 / srate;
        handles = set_Wd(handles,Wd);
    end
    if ~isempty(get(handles.Ws,'String'))
        Ws      = round(str2double(get(handles.Ws,'String')) * srate / 1000) * 1000 / srate;
        handles = set_Ws(handles,Ws);
    end
    if ~isempty(get(handles.Wp,'String'))
        Wp      = round(str2double(get(handles.Wp,'String')) * srate / 1000) * 1000 / srate;
        handles = set_Wp(handles,Wp);
    end
    if ~isempty(get(handles.Wv,'String'))
        Wv      = round(str2double(get(handles.Wv,'String')) * srate / 1000) * 1000 / srate;
        handles = set_Wv(handles,Wv);
    end
end
end

function handles = reset_parameters(handles)
handles = set_fname(handles,[]);
handles = set_nchannels(handles,[]);
handles = set_nframes(handles,[]);
handles = set_ntrials(handles,[]);
handles = set_start_time(handles,[]);
handles = set_srate(handles,[]);
handles = set_Wr(handles,[]);
handles = set_Wd(handles,[]);
handles = set_Ws(handles,[]);
handles = set_Wp(handles,[]);
handles = set_Wv(handles,[]);
handles = set_segindicator(handles,[]);
handles = copy_parameters(handles,true);
if isfield(handles,'EEG')
    handles = rmfield(handles,'EEG');
end
end

function handles = copy_parameters(handles,reset)
if (nargin<2) || ~reset
    T0                     = get(handles.start_time,'String');
    fs                     = get(handles.srate,'String');
    Wr                     = get(handles.Wr,'String');
    Wd                     = get(handles.Wd,'String');
    Ws                     = get(handles.Ws,'String');
    Wp                     = get(handles.Wp,'String');
    Wv                     = get(handles.Wv,'String');
    handles.EEG.start_time = str2double(T0);
    handles.EEG.srate      = str2double(fs);
    handles.EEG.Wr         = str2double(Wr);
    handles.EEG.Wd         = str2double(Wd);
    handles.EEG.Ws         = str2double(Ws);
    handles.EEG.Wp         = str2double(Wp);
    handles.EEG.Wv         = str2double(Wv);
    set(handles.T0_param,'String',T0);
    set(handles.fs_param,'String',fs);
    set(handles.Wr_param,'String',Wr);
    set(handles.Wd_param,'String',Wd);
    set(handles.Ws_param,'String',Ws);
    set(handles.Wp_param,'String',Wp);
    set(handles.Wv_param,'String',Wv);
else
    set(handles.T0_param,'String','');
    set(handles.fs_param,'String','');
    set(handles.Wr_param,'String','');
    set(handles.Wd_param,'String','');
    set(handles.Ws_param,'String','');
    set(handles.Wp_param,'String','');
    set(handles.Wv_param,'String','');
end
end

function handles = set_segindicator(handles,segindicator)
set(handles.segindicator,'String',string(segindicator))
end

function handles = set_segpnts(handles,segpnts)
if ~isempty(get(handles.ntrials,'String'))
    if iscell(segpnts) && (length(segpnts)==str2double(get(handles.ntrials,'String')))
        handles.EEG.segpnts = segpnts;
        handles             = set_segindicator(handles,true);
        handles             = copy_parameters(handles,false);
        handles             = set_time(handles);
    else
        handles = set_segindicator(handles,false);
        errordlg('Unrecognized segmentation structure','Loading Error')
    end
else
    handles = set_segindicator(handles,false);
    errordlg('#trials is missing!','Segmentation Error')
end
end

function handles = set_time(handles)
if ~isempty(get(handles.start_time,'String')) && ~isempty(get(handles.srate,'String')) ...
        && ~isempty(get(handles.nframes,'String'))
    start_time       = str2double(get(handles.start_time,'String'));
    srate            = str2double(get(handles.srate,'String'));
    nframes          = str2double(get(handles.nframes,'String'));
    time             = start_time + (0:nframes-1) * 1000 / srate;
    handles.EEG.time = time;
else
    errordlg('Some info entries are missing!','Time Error')
end
end
