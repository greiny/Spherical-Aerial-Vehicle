function varargout = droneAnim4(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @droneAnim4_OpeningFcn, ...
                   'gui_OutputFcn',  @droneAnim4_OutputFcn, ...
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


% --- Executes just before droneAnim4 is made visible.
function droneAnim4_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to droneAnim4 (see VARARGIN)

% Choose default command line output for droneAnim4
handles.output = hObject;

surf(handles.axes1);
rotate3d on;

surf(handles.axes2);
rotate3d on;

% Update handles structure
guidata(hObject, handles);
handles.j     = 1; % Initialize animation frame counting variable
handles.skipFlag = 0;
% Move GUI to center of screen
movegui('center')
% This sets up the initial plot - only do when we are invisible
% so window can get raised using droneAnim4.
A = evalin('base', 'yout');
% Generate the geometry used to draw the Sperical drone
r = 1; d = 1; h = 0.8; %inches: rotor dia., quad motor distance from 
% center of mass, and rotor height above arms (entirely cosmetic)
a = 0.8; b = 0.1; c = 0.8; i=0.3; g=0.2; % Top , bottom, and sides lengths of hub (?)

% Construct rotor representations
MOTOR = [0  0 h].';% m1 rotor center [X Y Z]
MOTORr = circlePoints(MOTOR, r, 10); MOTORr = [MOTORr MOTORr(:,1)]; % Rotor blade circles
% Motors connecting to center of blade circles
mMOTOR = [0,0;
      0,0;
      h,-b];
% Construct body plot points
bNS = [ 0, 0;
        c-g,  -c+g;
        g,  g]; %For drawing the body "X" shape
bEW = [ 0,  0;
        d, -d;
        i,  i];
% Body (HUB) Squares
Nfin = [ a+b, b, b; 0, 0, 0;  0,0,-c]; % By the way... north east is actually north west from above since x is north and y is east :P
Sfin = [ -a-b,-b, -b; 0, 0, 0; 0, 0, -c];
Efin = [0,0,0; -b,-a-b,-b; 0,0,-c ];
Wfin = [0,0,0; b,a+b,b; 0,0,-c ];

phi = A(1,4);
the = A(1,5);
psi = A(1,6);

% ROTATION MATRIX --- ZYX ROTATION (R = Rib)
R = [cos(psi)*cos(the) cos(psi)*sin(the)*sin(phi)-sin(psi)*cos(phi) cos(psi)*sin(the)*cos(phi)+sin(psi)*sin(phi);
       sin(psi)*cos(the) sin(psi)*sin(the)*sin(phi)+cos(psi)*cos(phi) sin(psi)*sin(the)*cos(phi)-cos(psi)*sin(phi);
       -sin(the)         cos(the)*sin(phi)                            cos(the)*cos(phi)];

% Rotate body frame velocity vector
U = A(:,7);
V = A(:,8);
W = A(:,9);
Vi = zeros(length(A),3);
MvMax= max(sqrt(U.^2+V.^2+W.^2)); 
Vb = 3/MvMax*[U(1), V(1), W(1)]'; % Scale velocity
Vi(1,:) = R*Vb; % Rotate velocity vector to inertial frame for plotting

% Support for X-configuration nifty trick (Requires that quadModel
% structure is in base workspace)

% Rotate body parts Via Initialized R
MOTORrR = R*MOTORr;
mMOTOR = R*mMOTOR;
bNSR = R*bNS;
bEWR = R*bEW;
NfinR = R*Nfin;
SfinR = R*Sfin;
EfinR = R*Efin;
WfinR = R*Wfin;


% Plot the box rotation and ang. velocity and inertial frame velocity
% vector
axes(handles.axes1)
plot3(bNSR(1,:),bNSR(2,:),bNSR(3,:),'k','LineWidth',10) % Body Arm
hold on
plot3(bEWR(1,:),bEWR(2,:),bEWR(3,:),'k','LineWidth',10) % Body Arm
plot3(mMOTOR(1,:),mMOTOR(2,:),mMOTOR(3,:),'k','LineWidth',10) % MotorS
plot3(MOTORrR(1,:),MOTORrR(2,:),MOTORrR(3,:),'g') % blades
grey = [0.5 0.5 0.5];
ne  = fill3(NfinR(1,:),NfinR(2,:),NfinR(3,:),'c'); alpha(ne,0.8); 
nw  = fill3(SfinR(1,:),SfinR(2,:),SfinR(3,:),grey); alpha(nw,0.8); 
sw  = fill3(EfinR(1,:),EfinR(2,:),EfinR(3,:),grey); alpha(sw,0.8); 
se  = fill3(WfinR(1,:),WfinR(2,:),WfinR(3,:),grey); alpha(se,0.8);
axis square
xlabel('X')
ylabel('Y')
zlabel('Z')
xlim([-3 3])
ylim([-3 3])
zlim([-3 3])
grid on


omi = zeros(length(A),3); % Initialize omega inertiaF points array (L(A)x3)
P = A(:,1); Q = A(:,2); Rw = A(:,3);
MombMax = max(sqrt(P.^2+Q.^2+Rw.^2)); % Calculate max magnitude of omb
omb = 3/MombMax*[P,Q,Rw].'; % Store and scale current omega bodyF
omi(1,:) = R*omb(:,1); % Rotate omegab to inertiaF store in omegai array
qp1 = quiver3(0,0,0,omi(1,1),omi(1,2),omi(1,3),'ro');
qp2 = quiver3(0,0,0,Vi(1,1),Vi(1,2),Vi(1,3),'k');
hold off
% legend([qp1 qp2],'Angular Velocity','I-Frame Velocity','Location',...
%      'SouthOutside')
%  pp1 = plot3(omi(1,1),omi(1,2),omi(1,3),'b.');
%  pp2 = plot3(Vi(1,1),Vi(1,2),Vi(1,3),'k.');
%  legend([qp1 qp2 pp1 pp2],'Angular Velocity','I-Frame Velocity',...
%      'Ang. Vel. Track','I-Frame Vel. Track','Location','SouthOutside')
% hold off

minX = min(A(:,10));
minY = min(A(:,11));
maxX = max(A(:,10));
maxY = max(A(:,11));
maxZ = max(A(:,12));

% Plot the three dimensional trajectory of the box
axes(handles.axes2)
X = A(1:1,10); Y = A(1:1,11); Z = A(1:1,12);
scatter3(X,Y,Z,36,'blue')
hold on
fill3([minX-1 maxX+1 maxX+1 minX-1],...
    [minY-1 minY-1 maxY+1 maxY+1],...
    [0 0 0 0],'g');
%alpha(0.7);
xlabel('X (m)')
ylabel('Y (m)')
zlabel('Z (m)')
xlim([minX-1 maxX+1])
ylim([minY-1 maxY+1])
zlim([-0.1 maxZ+1])
axis square
grid on
hold off
view(handles.axes1,35,25);
view(handles.axes2,35,25);

% phi_A = A(:,4);
% the_A = A(:,5);
% psi_A = A(:,6);

% Choose default command line output for animationGUI2
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes droneAnim4 wait for user response (see UIRESUME)
% uiwait(handles.figure1);



% --- Outputs from this function are returned to the command line.
function varargout = droneAnim4_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in startStop.
function startStop_Callback(hObject, eventdata, handles)
% hObject    handle to startStop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plotAnim(hObject, handles) % Pass handle to start/stop button to plot function

% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1


function points = circlePoints(center, radius, numberOfPoints)
% Helper function for plotting points
% Inspired by "circle()" from Peter Corke's MATLAB Robotics Toolbox
c = center.'; % [x;y] location of center
r = radius;
n = numberOfPoints;
% compute points around the circumference
th = (0:n-1)'/n*2*pi; % angles coresponding to each point
x = r*cos(th) + c(1); % x part
y = r*sin(th) + c(2); % y part
points = [x,y].';
    if length(c) > 2
        z = ones(size(x))*c(3); % z part
        points = [x, y, z].';
    end



function frameSkips_Callback(hObject, eventdata, handles)
% hObject    handle to frameSkips (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of frameSkips as text
%        str2double(get(hObject,'String')) returns contents of frameSkips as a double
handles.frameSkipVal = str2double(get(hObject,'String'));
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function frameSkips_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frameSkips (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.frameSkipVal = str2double(get(hObject,'String'));
guidata(hObject,handles);

function plotAnim(button, handles)
hObject = button;
mode = get(hObject,'String');
handles = guidata(gcbo);
j = handles.j;
if strcmp(mode,'Start')
    set(hObject,'String','Stop')
    guidata(hObject,handles)
    A = evalin('base', 'yout');
    tout = evalin('base', 'tout');
    frameSkipVal = str2double(get(handles.frameSkips,'String'))+1; % Size of steps to take for plotting animation (1 plots 
    % every frame, 2 about every other, etc.
    
    r = 1; d = 1; h = 0.8; %inches: rotor dia., quad motor distance from 
    % center of mass, and rotor height above arms (entirely cosmetic)
    a = 0.8; b = 0.1; c = 0.8; i=0.3; g=0.2; % Top , bottom, and sides lengths of hub (?)

    % Construct rotor representations
    MOTOR = [0  0 h].';% m1 rotor center [X Y Z]
    MOTORr = circlePoints(MOTOR, r, 10); MOTORr = [MOTORr MOTORr(:,1)]; % Rotor blade circles
    % Motors connecting to center of blade circles
    mMOTOR = [0,0;
              0,0;
              h,-b];
    % Construct body plot points
    bNS = [ 0, 0;
            c-g, -c+g;
            g,  g]; %For drawing the body "X" shape
    bEW = [ 0,  0;
            d, -d;
            i,  i];
    % Body (HUB) Squares
    Nfin = [ a+b, b, b; 0, 0, 0;  0,0,-c]; % By the way... north east is actually north west from above since x is north and y is east :P
    Sfin = [ -a-b,-b, -b; 0, 0, 0; 0, 0, -c];
    Efin = [0,0,0; -b,-a-b,-b; 0,0,-c ];
    Wfin = [0,0,0; b,a+b,b; 0,0,-c ];

    phi = A(1,4);
    the = A(1,5);
    psi = A(1,6);

    % Rib
    R = [cos(psi)*cos(the) cos(psi)*sin(the)*sin(phi)-sin(psi)*cos(phi) cos(psi)*sin(the)*cos(phi)+sin(psi)*sin(phi);
           sin(psi)*cos(the) sin(psi)*sin(the)*sin(phi)+cos(psi)*cos(phi) sin(psi)*sin(the)*cos(phi)-cos(psi)*sin(phi);
           -sin(the)         cos(the)*sin(phi)                            cos(the)*cos(phi)];
    
    % Bring in and rotate body frame velocity vector
    U = A(:,7);
    V = A(:,8);
    W = A(:,9);
    Vi = zeros(length(A),3);
    MvMax= max(sqrt(U.^2+V.^2+W.^2)); 
    Vb = 3/MvMax*[U, V, W].'; % Scale velocity vector points for plotting
    Vi(1,:) = R*Vb(:,1); % Rotate velocity vector to inertial frame for plotting
    % Angular Velocity Vector
    P = A(:,1); 
    Q = A(:,2); 
    Rw = A(:,3);
    omi = zeros(length(A),3); % Initialize omega inertiaF points array (L(A)x3)
    MombMax = max(sqrt(P.^2+Q.^2+Rw.^2)); % Calculate max magnitude of omb
    omb = 3/MombMax*[P,Q,Rw].'; % Store and scale current omega bodyF
    omi(1,:) = R*omb(:,1); % Rotate omegab to inertiaF store in omegai array

    phi = A(:,4);
    the= A(:,5);
    psi = A(:,6);

    % Next we run through the points in the vector A as an
    % animation until we get a command to stop
    minX = min(A(:,10)); % meters to feet
    minY = min(A(:,11));
    maxX = max(A(:,10));
    maxY = max(A(:,11));
    maxZ = max(A(:,12));
    colors = jet(length(A(:,10))); % color the path for time info
    
    % Precompute R
    % ROTATION MATRIX BELOW --- ZYX ROTATION
    R = cell(length(A),1); % create empty cell array for R matrices
    for i = 1:length(A)
    R{i,1} = [cos(psi(i))*cos(the(i)) cos(psi(i))*sin(the(i))*sin(phi(i))-sin(psi(i))*cos(phi(i)) cos(psi(i))*sin(the(i))*cos(phi(i))+sin(psi(i))*sin(phi(i));
              sin(psi(i))*cos(the(i)) sin(psi(i))*sin(the(i))*sin(phi(i))+cos(psi(i))*cos(phi(i)) sin(psi(i))*sin(the(i))*cos(phi(i))-cos(psi(i))*sin(phi(i));
              -sin(the(i))         cos(the(i))*sin(phi(i))                            cos(the(i))*cos(phi(i))];
    end
    % set up a j vector to straighten out skipping issues
%     if (skipFlag == 1)
%         jSkip = 1:frameSkipVal:j;
%     end
   
    view(handles.axes2,35,25);
    while ( strcmp(get(hObject,'String'),'Stop'))
        handles = guidata(gcbo);
        j = handles.j;
        
        guidata(hObject,handles); % Update handles data
        Vi(j,:) = R{j,1}*Vb(:,j);

        MOTORrR = R{j,1}*MOTORr;
        bNSR = R{j,1}*bNS;
        bEWR = R{j,1}*bEW;
        % Rotate body parts Via Initialized R
        mMOTORr = R{j,1}*mMOTOR;
        NfinR = R{j,1}*Nfin;
        SfinR = R{j,1}*Sfin;
        EfinR = R{j,1}*Efin;
        WfinR = R{j,1}*Wfin;

        % Plot the quad rotation and ang. velocity and inertial frame velocity vector
        axes(handles.axes1)
        plot3(bNSR(1,:),bNSR(2,:),bNSR(3,:),'k','LineWidth',10)
        hold on
        plot3(bEWR(1,:),bEWR(2,:),bEWR(3,:),'k','LineWidth',10)
        plot3(MOTORrR(1,:),MOTORrR(2,:),MOTORrR(3,:),'g')
        plot3(mMOTORr(1,:),mMOTORr(2,:),mMOTORr(3,:),'k','LineWidth',10) % Motor
        grey = [0.5 0.5 0.5];
        ne  = fill3(NfinR(1,:),NfinR(2,:),NfinR(3,:),'c'); alpha(ne,0.8); % North East surface
        nw  = fill3(SfinR(1,:),SfinR(2,:),SfinR(3,:),grey); alpha(nw,0.8); % North West surface
        sw  = fill3(EfinR(1,:),EfinR(2,:),EfinR(3,:),grey); alpha(sw,0.8); % South West surface
        se  = fill3(WfinR(1,:),WfinR(2,:),WfinR(3,:),grey); alpha(se,0.8); % South East surface
        xlabel('X')
        ylabel('Y')
        zlabel('Z')
        xlim([-3 3])
        ylim([-3 3])
        zlim([-3 3])
        view(handles.axes1,35,25);
        % set(gca,'Color',[215/255 244/255 247/255]) % This is a technique
        % that can be used for setting the plot background color
        
%         P = A(j,1); Q = A(j,2); Rw = A(j,3);
%         M = sqrt(P^2+Q^2+Rw^2);
%         omb = 3/M*[P,Q,Rw].';%/M; % Scaling of angular velocity vector
        omi(j,:) = R{j,1}*omb(:,j);
        qp1 = quiver3(0,0,0,omi(j,1),omi(j,2),omi(j,3),'r');
        qp2 = quiver3(0,0,0,Vi(j,1),Vi(j,2),Vi(j,3),'k');
        axis square
        grid on
        hold off
        drawnow
       
        % Plot the three dimensional trajectory of the box
        if (j==1)
            cla(handles.axes2)
        else if (handles.skipFlag==1)
                cla(handles.axes2)
                X = A(1:frameSkipVal:j,10); 
                Y = A(1:frameSkipVal:j,11); 
                Z = A(1:frameSkipVal:j,12);
                axes(handles.axes2)
                hold on
                scatter3(X,Y,Z,36,colors(1:frameSkipVal:j,:)); 
            end
        end
        % Ordinary plotting sequence (one frame at a time using hold on to
        % create persistance
        X = A(j,10); Y = A(j,11); Z = A(j,12);
        axes(handles.axes2)
        hold on
        scatter3(X,Y,Z,36,colors(j,:));
        if (j == 1 || handles.skipFlag==1)
            fill3([minX-1 maxX+1 maxX+1 minX-1],...
                  [minY-1 minY-1 maxY+1 maxY+1],...
                  [0 0 0 0],'g'); % make a plane to represent the ground (Z = 0)
            alpha(0.7); % Makes the ground "see-through"
            handles.skipFlag = 0;
            guidata(hObject,handles)
        end
        
        xlabel('X (m)')
        ylabel('Y (m)')
        zlabel('Z (m)')
        xlim([minX-1 maxX+1])
        ylim([minY-1 maxY+1])
        zlim([-0.1 maxZ+1]) % Keep ground within view
        % These calls are necessarry to avoid stutter in the frames if adjusting slider while animation is running
        % Set the current vie
        axis square
%         axis equal
%         xlimz = xlim;
%         ylimz = ylim;
%                     fill3([xlimz(1) xlimz(2) xlimz(2) xlimz(1)],...
%                   [ylimz(1) ylimz(1) ylimz(2) ylimz(2)],...
%                   [0 0 0 0],'k'); % make a plane to represent the ground (Z = 0)
%             alpha(0.9); % Makes the ground "see-through"
        grid on
        drawnow
        hold off
%         if (firstRunFlag==1)
%             firstRunFlag = 0;
%         end
        if (j == size(A,1)) % If we've reached the last group of data
            handles.j = 1; % Reset j to 1
            set(hObject,'String','Start'); % Reset button state
            guidata(hObject,handles);
            return % exit the function
            
        else if (j+frameSkipVal<size(A,1)) % If the next index won't exceed total frames...
        j = j+frameSkipVal; % increase the index by frame skip value
            else
            j = size(A,1); % We're almost out of data so to the last value for the final plot
            end
        end
        set(handles.simTime,'String',num2str(tout(j))); % update time display
        handles = guidata(gcbo); % update handles structure again...
        if (handles.skipFlag == 1)
            guidata(hObject,handles) % this nightmare is used to catch skip events which otherwise tend to go unheeded.
            return
        else
        handles.j = j; % store j as a handles value for persistance
        guidata(hObject,handles); 
        end
    end
    else if strcmp(get(hObject,'String'),'Stop')
        set(hObject,'String','Start');
        end
end


% --- Executes during object creation, after setting all properties.
function simTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to simTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function axes2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes2


% --- Executes on button press in XYview.
function XYview_Callback(hObject, eventdata, handles)
% hObject    handle to XYview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
view(handles.axes1,0,90);
view(handles.axes2,0,90);
drawnow
guidata(hObject,handles);


% --- Executes on button press in XZview.
function XZview_Callback(hObject, eventdata, handles)
% hObject    handle to XZview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
view(handles.axes1,0,0);
view(handles.axes2,0,0);
drawnow
guidata(hObject,handles);


% --- Executes on button press in YZview.
function YZview_Callback(hObject, eventdata, handles)
% hObject    handle to YZview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
view(handles.axes1,90,0);
view(handles.axes2,90,0);
drawnow
guidata(hObject,handles);


% --- Executes on button press in defaultView.
function defaultView_Callback(hObject, eventdata, handles)
% hObject    handle to defaultView (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
view(handles.axes1,35,25);
view(handles.axes2,35,25);
drawnow
guidata(hObject,handles);


% --- Executes on button press in timeSkipButton.
function timeSkipButton_Callback(hObject, eventdata, handles)
% hObject    handle to timeSkipButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.startStop,'String','Start');
tout = evalin('base', 'tout');
TimeSkip = str2double(get(handles.timeSkipEditBox,'String'));
handles.j = find(tout<=TimeSkip, 1,'last');
handles.skipFlag = 1;
% if handles.doubleRunFlag == 0
%     timeSkipButton_Callback(hObject, eventdata, handles) % run it again
% end
% handles.doubleRunFlag = 1;
guidata(hObject,handles);
% yes i'm seriously resorting to this...





function timeSkipEditBox_Callback(hObject, eventdata, handles)
% hObject    handle to timeSkipEditBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of timeSkipEditBox as text
%        str2double(get(hObject,'String')) returns contents of timeSkipEditBox as a double
handles.TimeSkip = str2double(get(hObject,'String'));
handles.skipFlag = 0;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function timeSkipEditBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to timeSkipEditBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.TimeSkip = 0;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
