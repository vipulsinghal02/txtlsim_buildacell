function varargout = txtl_plot_gui(varargin)
% TXTL_PLOT_GUI M-file for txtl_plot_gui.fig
%      TXTL_PLOT_GUI, by itself, creates a new TXTL_PLOT_GUI or raises the existing
%      singleton*.
%
%      H = TXTL_PLOT_GUI returns the handle to a new TXTL_PLOT_GUI or the handle to
%      the existing singleton*.
%
%      TXTL_PLOT_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TXTL_PLOT_GUI.M with the given input arguments.
%
%      TXTL_PLOT_GUI('Property','Value',...) creates a new TXTL_PLOT_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before txtl_plot_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to txtl_plot_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help txtl_plot_gui

% Last Modified by GUIDE v2.5 06-Jan-2013 21:25:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @txtl_plot_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @txtl_plot_gui_OutputFcn, ...
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


% --- Executes just before txtl_plot_gui is made visible.
function txtl_plot_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to txtl_plot_gui (see VARARGIN)

% Choose default command line output for txtl_plot_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% get default data structrure for plot
defaultdataGroups = txtl_getDefaultPlotDataStruct();

switch nargin
    
    case 5
        simData = varargin{1}; 
        modelObj = varargin{2};
        t_ode = simData.Time;
        x_ode = simData.Data;
        dataGroups = defaultdataGroups;
    case 6
        t_ode = varargin{1};
        x_ode = varargin{2};
        modelObj = varargin{3};
        dataGroups = defaultdataGroups;
    case 7
        t_ode = varargin{1};
        x_ode = varargin{2};
        modelObj = varargin{3};
        dataGroups = varargin{4};
    otherwise
        error('');
end


processedData = txtl_plot(t_ode,x_ode,modelObj,dataGroups,handles);
% set up the checkboxes
numberOfProteins = size(processedData{1}{1},2)-1; % first column is time
numOfcheckBoxes = 10;

for k=1:numOfcheckBoxes
    if k <= numberOfProteins
        eval(sprintf('set(handles.checkbox%d,''String'',processedData{1}{1}{k+1})',k));
        eval(sprintf('set(handles.checkbox%d,''Value'',1)',k));
        eval(sprintf('set(handles.checkbox%d,''Visible'',''on'')',k));
    else
        eval(sprintf('set(handles.checkbox%d,''Visible'',''off'')',k));
    end 
end

% UIWAIT makes txtl_plot_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = txtl_plot_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1
setLineVisibility(1,get(hObject,'Value'),handles) 


% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2
setLineVisibility(2,get(hObject,'Value'),handles) 


% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox3
setLineVisibility(3,get(hObject,'Value'),handles) 


% --- Executes on button press in checkbox4.
function checkbox4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox4
setLineVisibility(4,get(hObject,'Value'),handles) 


% --- Executes on button press in checkbox5.
function checkbox5_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox5
setLineVisibility(5,get(hObject,'Value'),handles) 

% --- Executes on button press in checkbox6.
function checkbox6_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox6
setLineVisibility(6,get(hObject,'Value'),handles) 

% --- Executes on button press in checkbox7.
function checkbox7_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox7
setLineVisibility(7,get(hObject,'Value'),handles) 

% --- Executes on button press in checkbox8.
function checkbox8_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox8
setLineVisibility(8,get(hObject,'Value'),handles) 

% --- Executes on button press in checkbox9.
function checkbox9_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox9
setLineVisibility(9,get(hObject,'Value'),handles) 

% --- Executes on button press in checkbox10.
function checkbox10_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox10
setLineVisibility(10,get(hObject,'Value'),handles) 

function setLineVisibility(numOfLine,status,handles)
    lines = get(handles.genePlot,'Children');
    lines = flipud(lines);
    if status 
        set(lines(numOfLine),'Visible','on');
    else
        set(lines(numOfLine),'Visible','off');
    end
    
