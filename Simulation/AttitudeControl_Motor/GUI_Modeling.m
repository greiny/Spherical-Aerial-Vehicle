function varargout = GUI_Modeling(varargin)

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_Modeling_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_Modeling_OutputFcn, ...
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

function GUI_Modeling_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_Modeling (see VARARGIN)

imshow('Drone_Graphic.jpg'); % Displays default graphic before user selects one 

movegui('center') % Moves the current GUI to the center of screen

% Choose default command line output for GUI_Modeling
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_Modeling wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = GUI_Modeling_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function drone_m_edit_Callback(hObject, eventdata, handles)
% hObject    handle to drone_m_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of drone_m_edit as text
%        str2double(get(hObject,'String')) returns contents of drone_m_edit as a double


% --- Executes during object creation, after setting all properties.
function drone_m_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to drone_m_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function drone_ix_edit_Callback(hObject, eventdata, handles)
% hObject    handle to drone_ix_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of drone_ix_edit as text
%        str2double(get(hObject,'String')) returns contents of drone_ix_edit as a double


% --- Executes during object creation, after setting all properties.
function drone_ix_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to drone_ix_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function drone_iy_edit_Callback(hObject, eventdata, handles)
% hObject    handle to drone_iy_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of drone_iy_edit as text
%        str2double(get(hObject,'String')) returns contents of drone_iy_edit as a double


% --- Executes during object creation, after setting all properties.
function drone_iy_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to drone_iy_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function drone_iz_edit_Callback(hObject, eventdata, handles)
% hObject    handle to drone_iz_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of drone_iz_edit as text
%        str2double(get(hObject,'String')) returns contents of drone_iz_edit as a double


% --- Executes during object creation, after setting all properties.
function drone_iz_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to drone_iz_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function jx_Callback(hObject, eventdata, handles)
% hObject    handle to jx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of jx as text
%        str2double(get(hObject,'String')) returns contents of jx as a double


% --- Executes during object creation, after setting all properties.
function jx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to jx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function jy_Callback(hObject, eventdata, handles)
% hObject    handle to jy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of jy as text
%        str2double(get(hObject,'String')) returns contents of jy as a double


% --- Executes during object creation, after setting all properties.
function jy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to jy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function jz_Callback(hObject, eventdata, handles)
% hObject    handle to jz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of jz as text
%        str2double(get(hObject,'String')) returns contents of jz as a double


% --- Executes during object creation, after setting all properties.
function jz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to jz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in save.
function save_Callback(hObject, eventdata, handles)
% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get the handles of variables desired
% Converts the string representations of the variables into numeric values


% Calculate and convert everything here to go from g/cm --> kg/m
    
motor_m= str2num(get(handles.drone_m_edit, 'String'));
motor_ix = str2num(get(handles.drone_ix_edit, 'String'));
motor_iy = str2num(get(handles.drone_iy_edit, 'String'));
motor_iz = str2num(get(handles.drone_iz_edit, 'String'));

h1 = str2num(get(handles.h1_edit, 'String'));
h2 = str2num(get(handles.h2_edit, 'String'));
h3 = str2num(get(handles.h3_edit, 'String'));
D = str2num(get(handles.D_edit, 'String'));
FA = str2num(get(handles.FA_edit, 'String'));

cl = str2num(get(handles.cl_edit, 'String'));
cq = str2num(get(handles.cq_edit, 'String'));
Alpha = str2num(get(handles.Alpha_edit, 'String'));
Beta = str2num(get(handles.Beta_edit, 'String'));

b = str2num(get(handles.b_edit, 'String'));
E = str2num(get(handles.E_edit, 'String'));
T = str2num(get(handles.tConstant, 'String')); % Seconds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%값 저장하기

g = 9.81; % m/s^2
mass = motor_m;
Ix = motor_ix;
Iy = motor_iy;
Iz = motor_iz;
Ib = [Ix 0 0; 0 Iy 0; 0 0 Iz];
Ibinv = [1/Ix 0 0; 0 1/Iy 0; 0 0 1/Iz];
% dhh matrix for "Plus" Configuration
plusConfig = 1;

% Saves all the variables from above into a structure called "quadModel"
droneModel = struct('g',(g),'mass',(mass),'Ib',(Ib),'Ibinv',(Ibinv),...
    'Ix',(Ix),'Iy',(Iy),'Iz',(Iz)','h1',(h1)','h2',(h2)','h3',(h3),...
    'D',(D),'E',(E)','T',(T),'Alpha',(Alpha),'Beta',(Beta),...
    'FA',(FA),'cq',(cq),'cl',(cl),'b',(b), 'plusConfig',(plusConfig));
   


% Assigns the structure to the active workspace
% assignin('WS','name',V) assigns the variable 'name' in the workspace WS the value V.
% assignin('base','quadModel',quadModel);

% Saves the structure "quadModel" into the working directory
uisave('droneModel','droneModel_+');

guidata(hObject, handles);

% --- Executes on button press in clear_button.
function clear_button_Callback(hObject, eventdata, handles)
% hObject    handle to clear_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Clear all the text boxes
set(handles.drone_m_edit,'String','');
set(handles.drone_ix_edit,'String','');
set(handles.drone_iy_edit,'String','');
set(handles.drone_iz_edit,'String','');

set(handles.h1_edit,'String','');
set(handles.h2_edit,'String','');
set(handles.h3_edit,'String','');
set(handles.D_edit,'String','');
set(handles.FA_edit,'String','');

set(handles.cl_edit,'String','');
set(handles.cq_edit,'String','');
set(handles.Alpha_edit,'String','');
set(handles.Beta_edit,'String','');

set(handles.b_edit,'String','');
set(handles.E_edit,'String','');
set(handles.tConstant,'String','');



% --- Executes on button press in loadButton.
function loadButton_Callback(hObject, eventdata, handles)
% hObject    handle to loadButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% This Load button is used when the user wishes to change only a few
% values of an existing structure (that is already saved). Instead of
% re-typing all the quad parameters back into the GUI again only to change
% 1 or 2 values, they can load the values back in quickly and make the necessary
% changes, and re-save as required.

uiload; % Brings up standard load window, a structure is then selected

if exist('droneModel')

set(handles.drone_m_edit,'String',droneModel.mass);
set(handles.drone_ix_edit,'String',droneModel.Ix);
set(handles.drone_iy_edit,'String',droneModel.Iy);
set(handles.drone_iz_edit,'String',droneModel.Iz); 

set(handles.h1_edit,'String',droneModel.h1); 
set(handles.h2_edit,'String',droneModel.h2); 
set(handles.h3_edit,'String',droneModel.h3); 
set(handles.D_edit,'String',droneModel.D); 
set(handles.FA_edit,'String',droneModel.FA);

set(handles.cl_edit,'String',droneModel.cl); 
set(handles.cq_edit,'String',droneModel.cq); 
set(handles.Alpha_edit,'String',droneModel.Alpha); 
set(handles.Beta_edit,'String',droneModel.Beta); 

set(handles.b_edit,'String',droneModel.b);
set(handles.E_edit,'String',droneModel.E); 
set(handles.tConstant,'String',droneModel.T);


else
end

function tConstant_Callback(hObject, eventdata, handles)
% hObject    handle to tConstant (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tConstant as text
%        str2double(get(hObject,'String')) returns contents of tConstant as a double


% --- Executes during object creation, after setting all properties.
function tConstant_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tConstant (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function h1_edit_Callback(hObject, eventdata, handles)
% hObject    handle to h1_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of h1_edit as text
%        str2double(get(hObject,'String')) returns contents of h1_edit as a double


% --- Executes during object creation, after setting all properties.
function h1_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to h1_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function h2_edit_Callback(hObject, eventdata, handles)
% hObject    handle to h2_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of h2_edit as text
%        str2double(get(hObject,'String')) returns contents of h2_edit as a double


% --- Executes during object creation, after setting all properties.
function h2_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to h2_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function h3_edit_Callback(hObject, eventdata, handles)
% hObject    handle to h3_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of h3_edit as text
%        str2double(get(hObject,'String')) returns contents of h3_edit as a double


% --- Executes during object creation, after setting all properties.
function h3_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to h3_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function D_edit_Callback(hObject, eventdata, handles)
% hObject    handle to D_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of D_edit as text
%        str2double(get(hObject,'String')) returns contents of D_edit as a double


% --- Executes during object creation, after setting all properties.
function D_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to D_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Alpha_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Alpha_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Alpha_edit as text
%        str2double(get(hObject,'String')) returns contents of Alpha_edit as a double


% --- Executes during object creation, after setting all properties.
function Alpha_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Alpha_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Beta_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Beta_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Beta_edit as text
%        str2double(get(hObject,'String')) returns contents of Beta_edit as a double


% --- Executes during object creation, after setting all properties.
function Beta_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Beta_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function E_edit_Callback(hObject, eventdata, handles)
% hObject    handle to E_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of E_edit as text
%        str2double(get(hObject,'String')) returns contents of E_edit as a double


% --- Executes during object creation, after setting all properties.
function E_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to E_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function FA_edit_Callback(hObject, eventdata, handles)
% hObject    handle to FA_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FA_edit as text
%        str2double(get(hObject,'String')) returns contents of FA_edit as a double


% --- Executes during object creation, after setting all properties.
function FA_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FA_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cq_edit_Callback(hObject, eventdata, handles)
% hObject    handle to cq_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cq_edit as text
%        str2double(get(hObject,'String')) returns contents of cq_edit as a double


% --- Executes during object creation, after setting all properties.
function cq_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cq_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cl_edit_Callback(hObject, eventdata, handles)
% hObject    handle to cl_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cl_edit as text
%        str2double(get(hObject,'String')) returns contents of cl_edit as a double


% --- Executes during object creation, after setting all properties.
function cl_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cl_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function b_edit_Callback(hObject, eventdata, handles)
% hObject    handle to b_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of b_edit as text
%        str2double(get(hObject,'String')) returns contents of b_edit as a double


% --- Executes during object creation, after setting all properties.
function b_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to b_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
