%% This demonstration script generates a timeseries in the workspace
%
X = [ 0, 1, 0,-1, 0, 1, 0, 0]; % meters
Y = [ 0, 0, 1, 0,-1, 0, 0, 0]; % meters
Z = [ 10, 10, 10, 10, 10, 10, 10, 10]; % meters
t = [ 0, 5,10,15,20,25,30, 35]; % seconds
Psi = [0,90*pi/180,0,180*pi/180,0,0,0,90*pi/180]; % radians
path.x = timeseries(X,t);
path.y = timeseries(Y,t);
path.z = timeseries(Z,t);
path.psi = timeseries(Psi,t);
clear X Y Z t Psi
uisave('path', 'Path_Diamond')
%clear path