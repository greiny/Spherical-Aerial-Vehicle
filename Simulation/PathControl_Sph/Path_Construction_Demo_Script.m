% X = [ 0, 1, 0,-1, 0, 1, 0, 0]; % meters
% Y = [ 0, 0, 1, 0,-1, 0, 0, 0]; % meters
% Z = [ 10, 10, 10, 10, 10, 10, 10, 10]; % meters
% t = [ 0, 5,10,15,20,25,30, 35]; % seconds
% %Psi = [0,90*pi/180,0,180*pi/180,0,0,0,90*pi/180]; % radians
% path.x = timeseries(X,t);
% path.y = timeseries(Y,t);
% path.z = timeseries(Z,t);
% %path.psi = timeseries(Psi,t);
% %clear X Y Z t Psi
% clear X Y Z t
% uisave('path', 'Path_Diamond')
% %clear path

% <Circular Path Construction Script>
% Direction is CW by default.
% The script will automatically assign time required for each linear
% path(s) based on user configurable variables below.

r = 3;          % Radius of the circle in meters
num_ver = 100;  % Number of vertex for approximation
t_tgt = 20;     % Target time to complete circuit in seconds
alt_tgt = 10;   % Target altitude in meters

X = 0:num_ver; %
Y = 0:num_ver; %
Z = 0:num_ver; % Constructing empty matrix for coordinates assignment
t = 0:(t_tgt/num_ver):t_tgt; % Allocation of time for each path

for i = 1:num_ver + 1
    X(i) = r * cos(2 * pi * i / num_ver) - r;
    Y(i) = r * sin(2 * pi * i / num_ver);
    Z(i) = alt_tgt;
    %Z(i) = 20 * (i / num_ver);
end

path.x = timeseries(X,t);
path.y = timeseries(Y,t);
path.z = timeseries(Z,t);
clear X Y Z t
uisave('path', 'Path_Circle_Approx')