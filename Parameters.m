%% GUI functions
function varargout = Parameters(varargin)
% PARAMETERS MATLAB code for Parameters.fig
%      PARAMETERS, by itself, creates a new PARAMETERS or raises the existing
%      singleton*.
%
%      H = PARAMETERS returns the handle to a new PARAMETERS or the handle to
%      the existing singleton*.
%
%      PARAMETERS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PARAMETERS.M with the given input arguments.
%
%      PARAMETERS('Property','Value',...) creates a new PARAMETERS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Parameters_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Parameters_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Parameters

% Last Modified by GUIDE v2.5 03-Oct-2018 03:14:35

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Parameters_OpeningFcn, ...
                   'gui_OutputFcn',  @Parameters_OutputFcn, ...
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


% --- Executes just before Parameters is made visible.
function Parameters_OpeningFcn(hObject, eventdata, handles, varargin)
% Choose default command line output for Parameters
handles.output = hObject;
get_actual_parameters(handles)

% Update handles structure
guidata(hObject, handles);




% --- Outputs from this function are returned to the command line.
function varargout = Parameters_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;





%% Help Functions
% DataBlock Help
function DataBlock_Help_Callback(hObject, eventdata, handles)
helpdlg(['DataBlock Size determines the  largest length of each channel ' ... 
    'of the raw data package to process at once. The value of this parameter ' ...
    'depends on user''s available memory for the program'],'DataBlock Parameter');

% Normalization Help
function Normalization_Help_Callback(hObject, eventdata, handles)
helpdlg(['This enables normalization of the amplitude of Raw Data. ' ... 
    'Also, it disables the signal scaling done by Scale parameter. ' ...
    'This option is useful to make a comparision of signals acquired with ' ...
    ' different or unknown amplifying values.'],'Normalization Parameter');

% Ps Help
function Ps_help_Callback(hObject, eventdata, handles)
helpdlg(['It determines the size of feature points on the Features Plot. ' ...
    'This option is useful to get a better observation of feature points ' ...
    'when they are many. '],'Point Size Parameter');

% BlockTime Help
function BlockTime_Help_Callback(hObject, eventdata, handles)
helpdlg(['It determines the maximum waiting time for plotting figures. While ' ... 
    'plotting, other functions are blocked to avoid errors. Although, it is possible '  ...
    'to disable this feature by setting 0 this parameter'],'Block Time Parameter');

% Density Resolution Help
function Density_Resolution_Help_Callback(hObject, eventdata, handles)
helpdlg(['The density plot is a 3D histogram of the Features points 2D space. ' ... 
    'This space is divided in cells forming a n*m grid. These 2 parameters are ' ...
    'the number of rows and columns of that grid. The more cells has the grid, the ' ...
    'more resolution is obtained, but less discrimination. It is suggested to set '...
    'rows and columns to the same value (square grid).'],'Density Resolution Parameter');

% Record Data Help
function Record_Data_Help_Callback(hObject, eventdata, handles)
helpdlg(['Destination Folder option is only available after loading a file. ' ... 
    'Also, the configuration of this parameter only works for the current file. ' ...
    'After loading a new one, it is necessary to specify again the directory or let ' ...
    'it remain as default.'],'Destination Folder');




%% Parameters  general functions

% Set parameters to actual user values (actual GUI parameter values) and refresh GUI values 
function Save_parameters_Callback(hObject, eventdata, handles)
set_user_parameters(handles);
get_actual_parameters(handles)

% Close GUI without setting parameters not previously saved.
function Cancel_parameters_Callback(hObject, eventdata, handles)
close(Parameters);

% Set all parameters to default and refresh GUI parameter valures.
function Default_parameter_Callback(hObject, eventdata, handles)
ovw=questdlg('Do you want to reset all parameters to default?' ,'Default Parameters','Yes','No','No');
if strcmp(ovw,'Yes')
    set_default_parameters
    get_actual_parameters(handles)
end

%Refresh Parameters Values to actual
function Refresh_Callback(hObject, eventdata, handles)
get_actual_parameters(handles)

% Change Record Directory Path
function Change_Folder_Callback(hObject, eventdata, handles)
SelPath=uigetdir;
SelPath=[SelPath '\'];
set(handles.RecordPathName_parameter,'String',SelPath);

%Set parameters values of GUI to actual parameters value
function get_actual_parameters(handles)
global handles2Spk
global Channels
global Fs 
global DataBlock 
global Scale
global SL   
global norm 
global GcaColor 
global WfsColor
global SelColor
global RawColor
global BlockTime 
global RecordPathName 
global Ls 
global ps 
global DensityRes 
global ISI

set(handles2Spk.Channels,'String',num2str(Channels));
set(handles.Channels_parameter,'String',num2str(Channels));
set(handles.Fs_parameter,'String',num2str(Fs));
set(handles.DataBlock_parameter,'String',num2str(DataBlock));
set(handles.Scale_parameter,'String',num2str(Scale));
set(handles.SL_parameter,'String',num2str(SL));
set(handles.norm_parameter,'Value',norm);
set(handles.Gca_R_parameter,'String',num2str(GcaColor(1)));
set(handles.Gca_G_parameter,'String',num2str(GcaColor(2)));
set(handles.Gca_B_parameter,'String',num2str(GcaColor(3)));
set(handles.Gco_R_parameter,'String',num2str(WfsColor(1)));
set(handles.Gco_G_parameter,'String',num2str(WfsColor(2)));
set(handles.Gco_B_parameter,'String',num2str(WfsColor(3)));
set(handles.Sel_R_parameter,'String',num2str(SelColor(1)));
set(handles.Sel_G_parameter,'String',num2str(SelColor(2)));
set(handles.Sel_B_parameter,'String',num2str(SelColor(3)));
set(handles.Raw_R_parameter,'String',num2str(RawColor(1)));
set(handles.Raw_G_parameter,'String',num2str(RawColor(2)));
set(handles.Raw_B_parameter,'String',num2str(RawColor(3)));
set(handles.BlockTime_parameter,'String',num2str(BlockTime));
if isempty(RecordPathName)
    set(handles.RecordPathName_parameter,'String',RecordPathName);
else
    set(handles.RecordPathName_parameter,'String',RecordPathName);
    set(handles.Change_Folder,'Visible','on');
end
set(handles.Ls_parameter,'String',num2str(Ls));
set(handles.ps_parameter,'String',num2str(ps));
set(handles.ISI_parameter,'String',num2str(ISI));
set(handles.Density_R_parameter,'String',num2str(DensityRes(1)));
set(handles.Density_C_parameter,'String',num2str(DensityRes(2)));

%Set parameters to default value
function set_default_parameters

global Channels
global Fs 
global DataBlock 
global Scale
global SL   
global norm 
global GcaColor 
global WfsColor
global SelColor
global RawColor
global BlockTime 
global RecordPathName 
global Ls 
global ps 
global DensityRes 
global ISI

Channels=25;
Fs=30030;
DataBlock=801000;
Scale=2.048;
SL=16;
norm=0;
GcaColor=[0 0.15 0.35];
WfsColor=[0.35 0.6 0.55];
SelColor=[1 1 1];
RawColor=[0.8 0.3 0.3];
BlockTime=7;
Ls=4;
ps=25;
DensityRes=[100 100];
ISI=2;
if isempty(RecordPathName)
    RecordPathName=[];
else
    set_default_RecordPath;
end

% Set parameters to GUI parameters values
function set_user_parameters(handles)
global Channels
global Fs 
global DataBlock 
global Scale
global SL   
global norm 
global GcaColor 
global WfsColor 
global SelColor
global RawColor
global BlockTime 
global RecordPathName 
global Ls 
global ps 
global DensityRes 
global ISI

Channels=str2double(get(handles.Channels_parameter,'String'));
Fs=str2double(get(handles.Fs_parameter,'String'));
DataBlock=str2double(get(handles.DataBlock_parameter,'String'));
Scale=str2double(get(handles.Scale_parameter,'String'));
SL=str2double(get(handles.SL_parameter,'String'));
norm=get(handles.norm_parameter,'Value');
GcaColor(1)=str2double(get(handles.Gca_R_parameter,'String'));
GcaColor(2)=str2double(get(handles.Gca_G_parameter,'String'));
GcaColor(3)=str2double(get(handles.Gca_B_parameter,'String'));
WfsColor(1)=str2double(get(handles.Gco_R_parameter,'String'));
WfsColor(2)=str2double(get(handles.Gco_G_parameter,'String'));
WfsColor(3)=str2double(get(handles.Gco_B_parameter,'String'));
SelColor(1)=str2double(get(handles.Sel_R_parameter,'String'));
SelColor(2)=str2double(get(handles.Sel_G_parameter,'String'));
SelColor(3)=str2double(get(handles.Sel_B_parameter,'String'));
RawColor(1)=str2double(get(handles.Raw_R_parameter,'String'));
RawColor(2)=str2double(get(handles.Raw_G_parameter,'String'));
RawColor(3)=str2double(get(handles.Raw_B_parameter,'String'));
BlockTime=str2double(get(handles.BlockTime_parameter,'String'));
RecordPathName=get(handles.RecordPathName_parameter,'String');
Ls=str2double(get(handles.Ls_parameter,'String'));
ps=str2double(get(handles.ps_parameter,'String'));
ISI=str2double(get(handles.ISI_parameter,'String'));
DensityRes(1)=str2double(get(handles.Density_R_parameter,'String'));
DensityRes(2)=str2double(get(handles.Density_C_parameter,'String'));

%Set Record Directory Path
function set_default_RecordPath
    global RecordPathName
    global fileType
    global PathName
    global FileName
    if fileType==3  %if file is recorded data
            RecordPathName=PathName;
        else
            RecordPathName=[PathName 'Saves_' FileName  '/']; 
    end







%% Objects Creation DO NOT EDIT
function Channels_parameter_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function Fs_parameter_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function DataBlock_parameter_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function Scale_parameter_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function SL_parameter_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function Ls_parameter_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ps_parameter_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit9_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit10_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit11_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function Gca_R_parameter_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function BlockTime_parameter_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit16_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function Gca_G_parameter_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function Gca_B_parameter_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function Gco_R_parameter_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function Gco_G_parameter_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function Gco_B_parameter_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function Density_R_parameter_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function Density_C_parameter_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit31_CreateFcn(hObject, eventdata, handles)
function ISI_parameter_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function Sel_R_parameter_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function Sel_G_parameter_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function Sel_B_parameter_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function Raw_R_parameter_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function Raw_G_parameter_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function Raw_B_parameter_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
