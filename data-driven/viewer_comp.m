function varargout = viewer_comp(varargin)
% VIEWER_COMP MATLAB code for viewer_comp.fig
%      VIEWER_COMP, by itself, creates a new VIEWER_COMP or raises the existing
%      singleton*.
%
%      H = VIEWER_COMP returns the handle to a new VIEWER_COMP or the handle to
%      the existing singleton*.
%
%      VIEWER_COMP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VIEWER_COMP.M with the given input arguments.
%
%      VIEWER_COMP('Property','Value',...) creates a new VIEWER_COMP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before viewer_comp_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to viewer_comp_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help viewer_comp

% Last Modified by GUIDE v2.5 23-Apr-2015 23:43:38

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @viewer_comp_OpeningFcn, ...
                   'gui_OutputFcn',  @viewer_comp_OutputFcn, ...
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


% --- Executes just before viewer_comp is made visible.
function viewer_comp_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to viewer_comp (see VARARGIN)

% Choose default command line output for viewer_comp
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes viewer_comp wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% wangyu starts to put data here for initialization

% handles.peaks = peaks(35);
% [x,y] = meshgrid(-8:0.5:8);
% r = sqrt(x.^2+y.*2);
% handles.current_data = handles.peaks;
% surf(handles.peaks);

assert(length(varargin)>=3);

handles.V = varargin{1};
handles.F = varargin{2};
handles.comps = varargin{3};

handles.A = 1;
handles.index_comps = 1;
handles.ncomps = size(handles.comps,3);

%[handles.t] = render_component(handles.V,handles.F,handles.dV); % do not
%call this function
%hold off;
%axis manual;

%handles.t = tsurf(handles.F,handles.V,'EdgeColor','none','FaceColor',[150,220,150]./255.*1.1,'FaceLighting','phong','SpecularStrength',1.);
[handles.t,~,~] = render_mesh(handles.V,handles.F);

% Choose default command line output for simple_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles); % this is important

% end of wangyu's code

% --- Outputs from this function are returned to the command line.
function varargout = viewer_comp_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1


% --- Executes during object creation, after setting all properties.
function uipanel1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function num_of_comps_Callback(hObject, eventdata, handles)
% hObject    handle to num_of_comps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of num_of_comps as text
%        str2double(get(hObject,'String')) returns contents of num_of_comps as a double

% wangyu
newNum = str2double(get(hObject,'String'));
%gatherAndUpdate(handles, newNum);
handles.index_comps = newNum;

set(handles.t,'Vertices',handles.V+handles.comps(:,:,handles.index_comps)*handles.A);
guidata(hObject, handles); % this is important


% --- Executes during object creation, after setting all properties.
function num_of_comps_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num_of_comps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function magnitude_Callback(hObject, eventdata, handles)
% hObject    handle to magnitude (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of magnitude as text
%        str2double(get(hObject,'String')) returns contents of magnitude as a double
% wangyu
handles.A = str2double(get(hObject,'String'));

set(handles.t,'Vertices',handles.V+handles.comps(:,:,handles.index_comps)*handles.A);
guidata(hObject, handles); % this is important

% --- Executes during object creation, after setting all properties.
function magnitude_CreateFcn(hObject, eventdata, handles)
% hObject    handle to magnitude (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end