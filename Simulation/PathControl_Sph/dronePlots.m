function [] = quadPlots(yout,tout)
close all; clc;
[~,column] = size(yout);
if column==24
    PC = true;
else PC = false;
end
close all; clc;
A = yout;
T = tout;
t = T;
P = A(:,1); % rad/s for now, converted before plotting
Q = A(:,2);
R = A(:,3);
Phi   = A(:,4); % In radians for now, converted before plotting
Theta = A(:,5);
Psi   = A(:,6);
Phi_cmd = A(:,19);
Theta_cmd = A(:,20);
Psi_cmd = A(:,21);
U = A(:,7);
V = A(:,8);
W = A(:,9);
X = A(:,10);
Y = A(:,11);
Z = A(:,12);
Z_cmd   = A(:,22);
mc1 = A(:,13) * 180/pi;
mc2 = A(:,14) * 180/pi;
mc3 = A(:,15) * 180/pi;
mc4 = A(:,16) * 180/pi;
mc5 = A(:,17);
mc6 = A(:,18);

if and((Z<=0),(W<=0))
    mc5 = 0;
    mc6 = 0;
end

if (PC==true)
X_cmd = A(:,23);
Y_cmd = A(:,24);
end

% Plots ___________________________________________________________________
figure
subplot(4,3,1)
plot(T,P,'b')
xlabel('Time (s)')
ylabel('Angular Velocity (deg/sec)')
xlim([min(t) max(t)])
title('P')
grid on

subplot(4,3,2)
plot(T,Q,'r')
xlabel('Time (s)')
ylabel('Angular Velocity (deg/sec)')
xlim([min(t) max(t)])
title('Q')
grid on

subplot(4,3,3)
plot(T,R,'g')
xlabel('Time (s)')
ylabel('Angular Velocity (deg/sec)')
xlim([min(t) max(t)])
title('R')
grid on

subplot(4,3,4)
plot(T,Phi*180/pi,'b')
hold on
plot(T,Phi_cmd*180/pi,'k--')
hold off
xlabel('Time (s)')
ylabel('Angle (deg)')
xlim([min(t) max(t)])
title('Phi')
grid on
 
subplot(4,3,5)
plot(T,Theta*180/pi,'r')
hold on
plot(T,Theta_cmd*180/pi,'k--')
hold off
xlabel('Time (s)')
ylabel('Angle(deg)')
xlim([min(t) max(t)])
title('Theta')
grid on

subplot(4,3,6)
plot(T,Psi*180/pi,'g')
hold on
plot(T,Psi_cmd*180/pi,'k--')
hold off
xlabel('Time (s)')
ylabel('Angle (deg)')
xlim([min(t) max(t)])
title('Psi')
grid on

subplot(4,3,7)
plot(T,U,'b')
xlabel('Time (s)')
ylabel('Velocity (m/s)')
xlim([min(t) max(t)])
title('U')
grid on

subplot(4,3,8)
plot(T,V,'r')
xlabel('Time (s)')
ylabel('Velocity (m/s)')
xlim([min(t) max(t)])
title('V')
grid on

subplot(4,3,9)
plot(T,W,'g')
xlabel('Time (s)')
ylabel('Velocity (m/s)')
xlim([min(t) max(t)])
title('W')
grid on

subplot(4,3,10)
plot(T,X,'b')
if (PC==true)
hold on
plot(T,X_cmd,'k--')
hold off
end
xlabel('Time (s)')
ylabel('Position (m)')
xlim([min(t) max(t)])
title('X')
grid on

subplot(4,3,11)
plot(T,Y,'r')
if (PC==true)
hold on
plot(T,Y_cmd,'k--')
hold off
end
xlabel('Time (s)')
ylabel('Position (m)')
xlim([min(t) max(t)])
title('Y')
grid on

subplot(4,3,12)
plot(T,Z,'g')
hold on
plot(T,Z_cmd,'k--')
hold off
xlabel('Time (s)')
ylabel('Position (m)')
xlim([min(t) max(t)])
title('Z')
grid on

% Create a new figure for fin angle
figure
subplot(6,1,1)
plot(T,mc1);
xlabel('Time (s)')
ylabel('Fin Command (degree)')
xlim([min(T) max(T)])
title('Fin 1')
grid on

subplot(6,1,2)
plot(T,mc2);
xlabel('Time (s)')
ylabel('Fin Command (degree)')
xlim([min(T) max(T)])
title('Fin 2')
grid on

subplot(6,1,3)
plot(T,mc3);
xlabel('Time (s)')
ylabel('Fin Command (degree)')
xlim([min(T) max(T)])
title('Fin 3')
grid on

subplot(6,1,4)
plot(T,mc4);
xlabel('Time (s)')
ylabel('Fin Command (degree)')
xlim([min(T) max(T)])
title('Fin 4')
grid on

subplot(6,1,5)
plot(T,mc5);
xlabel('Time (s)')
ylabel('Motor Throttle Command (%)')
xlim([min(T) max(T)])
title('Motor')
grid on

subplot(6,1,6)
plot(T,mc6);
xlabel('Time (s)')
ylabel('Motor Throttle Command (%)')
xlim([min(T) max(T)])
title('Motor')
grid on

end
