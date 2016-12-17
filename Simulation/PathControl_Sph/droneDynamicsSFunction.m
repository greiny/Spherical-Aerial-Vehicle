function droneDynamicsSFunction(block)
setup(block);
end
 
function setup(block)
 
  % Register the number of ports.
  block.NumInputPorts  = 7;
  block.NumOutputPorts = 12;
  
  % Set up the port properties to be inherited or dynamic.
  
  for i = 1:6;
  block.InputPort(i).Dimensions         = 1;
  block.InputPort(i).DirectFeedthrough = false;
  block.InputPort(i).SamplingMode      = 'Sample';
  end
  %------
  % This is the disturbance input
  block.InputPort(7).Dimensions        = 6; % torques x,y,z; forces x,y,z.
  block.InputPort(7).DirectFeedthrough = false;
  block.InputPort(7).SamplingMode      = 'Sample';
 
  %------
  for i = 1:12;
  block.OutputPort(i).Dimensions       = 1;
  block.OutputPort(i).SamplingMode     = 'Sample';
  end
 
  block.NumDialogPrms = 2;
  block.NumContStates = 12;
  
  block.SampleTimes = [0 0];
  block.SetAccelRunOnTLC(false);
  block.SimStateCompliance = 'DefaultSimState';
  block.RegBlockMethod('CheckParameters', @CheckPrms);
  block.RegBlockMethod('InitializeConditions', @InitializeConditions);
  block.RegBlockMethod('Outputs', @Outputs);
  block.RegBlockMethod('Derivatives', @Derivatives);
end
 
function CheckPrms(block)
     drone   = block.DialogPrm(1).Data;
     IC     = block.DialogPrm(2).Data;
end
 
function InitializeConditions(block)
% Initialize 7 States
 
IC = block.DialogPrm(2).Data;
 
% IC.Phi, IC.The, IC.Psi are in deg ... convert to rads
% IC.P, IC.Q, IC.R are in deg/s ... convert to rad/s
P = IC.P*pi/180; Q = IC.Q*pi/180; R = IC.R*pi/180; 
Phi = IC.Phi*pi/180; The = IC.The*pi/180; Psi = IC.Psi*pi/180;
U = IC.U; V = IC.V; W = IC.W; 
X = IC.X; Y = IC.Y; Z = IC.Z;
 
init = [P,Q,R,Phi,The,Psi,U,V,W,X,Y,Z];
 
for i=1:12
block.OutputPort(i).Data = init(i);
block.ContStates.Data(i) = init(i);
end
 
end
 
function Outputs(block)
    for i = 1:12;
        block.OutputPort(i).Data = block.ContStates.Data(i);
    end
end
 
function Derivatives(block)
drone = block.DialogPrm(1).Data;
 
% P Q R in units of rad/sec
P = block.ContStates.Data(1);
Q = block.ContStates.Data(2);
R = block.ContStates.Data(3);
% Phi The Psi in radians
Phi = block.ContStates.Data(4);
The = block.ContStates.Data(5);
Psi = block.ContStates.Data(6);
% U V W in units of m/s
U = block.ContStates.Data(7);
V = block.ContStates.Data(8);
W = block.ContStates.Data(9);
% X Y Z in units of m
X = block.ContStates.Data(10);
Y = block.ContStates.Data(11);
Z = block.ContStates.Data(12);
 
% r values in radians of control fin
r1 = block.InputPort(1).Data;
r2 = block.InputPort(2).Data;
r3 = block.InputPort(3).Data;
r4 = block.InputPort(4).Data;
Throttle_1 = block.InputPort(5).Data;
Throttle_2 = block.InputPort(6).Data;

Dist_tau = block.InputPort(7).Data(1:3);
Dist_F = block.InputPort(7).Data(4:6);
 
%------
% Calculate Thrust
rho = 1.25; % Air density
rr = 0.765; % Effective air flow
E = drone.E;
alpha = drone.Alpha;
beta = drone.Beta;
D = drone.D;
RPM_ratio = drone.b/100;

%----- 75percent of Max Thrust and 86~87 Throttle for hovering
RPM_1 = Throttle_1 * RPM_ratio; 
RPM_2 = Throttle_2 * RPM_ratio; 
Thrust_1 = E * (rr*(pi/2)*D^2*rho*alpha^2*(RPM_1*0.001)^(2*beta))^(1/3);
Thrust_2 = E * (rr*(pi/2)*D^2*rho*alpha^2*(RPM_2*0.001)^(2*beta))^(1/3);
Thrust = Thrust_1 + Thrust_2;

% Compute Rotation Matrix
% We use a Z-Y-X rotation
Rib = [cos(Psi)*cos(The) cos(Psi)*sin(The)*sin(Phi)-sin(Psi)*cos(Phi) cos(Psi)*sin(The)*cos(Phi)+sin(Psi)*sin(Phi);
       sin(Psi)*cos(The) sin(Psi)*sin(The)*sin(Phi)+cos(Psi)*cos(Phi) sin(Psi)*sin(The)*cos(Phi)-cos(Psi)*sin(Phi);
       -sin(The)         cos(The)*sin(Phi)                            cos(The)*cos(Phi)];
Rbi = Rib';
 
 
%--------
% Calculate moment in body axis
A = drone.FA; % Area of fin [m^2]
Ap = 2*pi*D^2/4; % Area of prop disc [m^2]
Cl = drone.cl; % C_L/AOA of Fin(rad)
Vair = sqrt(Thrust/(rho*Ap)) ; % flow velocity through Control Volume
F_L = 0.5*rho*Cl*A*Vair^2; % [N/rad]

Cq = drone.cq;
Torque_1 = Cq*rho*(RPM_1/60)^2*D^5; % [Nm]
Torque_2 = Cq*rho*(RPM_2/60)^2*D^5;
Torque = Torque_1-Torque_2;
Motor_q = [0;0;Torque;];

ge = [0; 0; -drone.g];
gb = Rbi*ge;

Dist_Fb = Rbi*Dist_F;

% Fin 1,2,3,4, CG
% Moment arm matrix from geographic center to each center of force
hx = [0, drone.h2, 0, -drone.h2, 0];
hy = [drone.h2, 0, -drone.h2, 0, 0];
hz = [-drone.h1, -drone.h1, -drone.h1, -drone.h1, -drone.h3];
fb(:,1) = [F_L*r1*cos(r1); 0; -F_L*r1*sin(r1)];
fb(:,2) = [0; -F_L*r2*cos(r2); -F_L*r2*sin(r2)];
fb(:,3) = [F_L*r3*cos(r3); 0; -F_L*r3*sin(r3)];
fb(:,4) = [0; -F_L*r4*cos(r4); -F_L*r4*sin(r4)];  
fb(:,5) = (drone.mass-0.2)*gb;

Mb = [0;0;0];
hb = zeros(3,3,5);

for i=1:5 
    hb(:,:,i) = [0, -hz(i), hy(i);
                 hz(i), 0, -hx(i);
                -hy(i), hx(i), 0];
    M(:,i) = hb(:,:,i) * fb(:,i);
    Mb = Mb + M(:,i);
end
Mb = Mb + Dist_tau + Motor_q;
   
drone.plusConfig = 0;

% Obtain dP dQ dR
Fb = [ F_L*(r1*cos(r1)) + F_L*(r3*cos(r3));
      -F_L*(r2*cos(r2)) - F_L*(r4*cos(r4));
      Thrust - F_L*(r1*sin(r1)) - F_L*(r2*sin(r2)) - F_L*(r3*sin(r3)) - F_L*(r4*sin(r4))];
  
omb_bi = [P; Q; R];
OMb_bi = [ 0,-R, Q;
           R, 0,-P;
          -Q, P, 0];
 
b_omdotb_bi = drone.Ibinv*(Mb-OMb_bi*drone.Ib*omb_bi);
H_Phi = [1,tan(The)*sin(Phi), tan(The)*cos(Phi);
         0,         cos(Phi),         -sin(Phi);
         0,sin(Phi)/cos(The),cos(Phi)/cos(The)];   
Phidot = H_Phi*omb_bi;
 
% Compute Velocity and Position derivatives of body frame
 
 
vb = [U;V;W];
b_dv = (1/drone.mass)*Fb+gb+Dist_Fb-OMb_bi*vb; % Acceleration in body frame (FOR VELOCITY)\
i_dp = Rib*vb; % Units OK SI: Velocity of body frame w.r.t inertia frame (FOR POSITION)
 
dP = b_omdotb_bi(1);
dQ = b_omdotb_bi(2);
dR = b_omdotb_bi(3);
dPhi = Phidot(1);
dTheta = Phidot(2);
dPsi = Phidot(3);
dU = b_dv(1);
dV = b_dv(2);
dW = b_dv(3);
dX = i_dp(1);
dY = i_dp(2);
dZ = i_dp(3);
% Rough rule to impose a "ground" boundary...could easily be improved...
if ((Z<=0) && (dZ<=0)) % better  version then before?
    dZ = 0;
    W = 0;
    Throttle_1 = 0;
    Throttle_2 = 0;
    block.ContStates.Data(12) = 0;
end
f = [dP dQ dR dPhi dTheta dPsi dU dV dW dX dY dZ].';
  %This is the state derivative vector
block.Derivatives.Data = f;
end


