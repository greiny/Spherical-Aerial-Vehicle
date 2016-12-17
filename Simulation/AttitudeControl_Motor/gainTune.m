S1 = stepinfo(yout(:,4))
S2 = stepinfo(yout(:,5))
S3 = stepinfo(yout(:,6))

s1t = S1.SettlingTime
s1t2 = ceil(s1t)
s2t = S2.SettlingTime
s2t2 = ceil(s2t)
s3t = S3.SettlingTime
s3t2 = ceil(s3t)

if s1t2 == 0
    t1 =0
else
    t1 = tout(s1t2)
end

if s2t2 == 0
    t2 =0
else
    t2 = tout(s2t2)
end

if s3t2 == 0
    t3 =0
else
    t3 = tout(s3t2)
end

o1 = S1.Overshoot
o2 = S2.Overshoot
o3 = S3.Overshoot

a = 5
b = 7


if t3<a && o3<b
    Psi_cmd = Psi_cmd + 10*pi/180
end

%if((t1<a) && (t2<a) && (t3<a))
%    if((o1<b) && (o2<b) && (o3<b))
%        if mod(N,3)==0
%            Phi_cmd = Phi_cmd + 10*pi/180
%        end
%        if mod(N,3)==1
%            The_cmd = The_cmd + 10*pi/180
%        end
%        if mod(N,3)==2
%            Psi_cmd = Psi_cmd + 10*pi/180
%        end
%        N = N+1
%    end
%end

if ((t1>=a) || (t2>=a))
    KP1 = KP1 + 0.1
    KD1 = KD1 - 0.05
end
if (((o1>=b) && (o1 < 1000)) || ((o2>=b) && (o2 < 1000)))
    KD1 = KD1 + 0.1
    KP1 = KP1 - 0.05
end
if t3>=a
    KP2 = KP2 + 0.1
    KD2 = KD2 - 0.05
end
if ((o3>=b) && (o3 < 1000))
    KD2 = KD2 + 0.1
    KP2 = KP2 - 0.05
end
total = [t1,t2,t3; o1,o2,o3;
        Phi_cmd*180/pi,The_cmd*180/pi,Psi_cmd*180/pi]

