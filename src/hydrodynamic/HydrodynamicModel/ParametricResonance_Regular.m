%--------------------------------------------------------------------------
%%% Constants %%%
%----------------
% Density
rho = 1000;     %[kg/m^3]
% Gravity
g=9.81;         %[m/s^2]

%--------------------------------------------------------------------------
%%% Geometry  %%%
%----------------
%Diameter
D=37.2;         %[m]
%Area
A=pi/4*D^2;     %[m^2]
%Draft
Df=198.1;       %[m]
%Metacentric height
GM=10.1;        %[m]
%Displaced water volume
DeltaV=A*Df;
%Distance from bottom to CoG
KG=89;
%Distance of Centre of gravity from mean free surface
Hg=Df-KG;
%
l=212.9;
% Radius of gyration
kyy=29.2;       %[m]


%-----------------------------------------------------------------------
%Initial Positions
Heave1=0;
Pitch1=0;

%--------------------------------------------------------------------------
%%% Inertia %%%
%--------------
m = 215872*10^3;      %Msss[kg]
I55=m/32*(2*D^2+3*l^2+32*kyy^2);  %Momment of Inertia
A33 = 1/12*pi*rho*D^3;   % Added mass heave[kg]
A55 = 1/12*pi*rho*D^2*(KG^3+(Df-KG)^3); % Added mass pitch[kgm]
I = [m+A33,I55+A55]';  %Inertia vector for Simulink

%------------------------------------------------------------------------
%%% Damping %%%%
%---------------
w3=sqrt(rho*g*A/(m+A33));  %Heave natural frequency
w5=sqrt(rho*g*DeltaV*GM/(I55+A55)); %Pitch natural frequency
e3=0.012;  %Damping ratio heave (from Jingrui et al )
e5=0.019;   %Damping ratio pitch (from Jingrui et al )
C3=2*w3*e3*(m+A33);  %Heave damping
C5=2*w5*e5*(I55+A55);  %Pitch damping
C=[C3,C5]';  %Damping vector for Simulink


%--------------------------------------------------------------------------
%Time
T=1/(w5/(2*pi));
dt =T/100;
EndTime=floor(18000/dt)*dt;
Time=0:dt:EndTime;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%--------------------------------------------------------------------------
%%% Input Waves %%%%

Amplitude = 3.3; % Wave amplitude [m]    
F_PR=w5*2/(2*pi); %Parametric resonance 
percent = 0.9; %Percentage of parametric frequency
Frequency = F_PR*percent; % Wave frequency [Hz]
Period=1/Frequency;
w= 2*pi*Frequency; %[rad/s]

wave = Amplitude*sin(w*Time); % Create input wave signal[m]


%Ramp to ogradually increase wave amplitude to reduce transient effects
Ramp = ones(1,length(Time));
for ri=1:floor(10/dt)
    Ramp(ri)=0;
end

for ri = floor(10/dt)+1:floor(300/dt)
Ramp(ri)=Ramp(ri)-((floor(300/dt)-(ri-floor(10/dt)))/(floor(300/dt)));
end
wave=wave.*Ramp;
Wave=[Time',wave'];  %Wave vector for Simulink

%--------------------------------------------------------------------------
%%% Excitation force %%%
%-----------------------
%Load Excitation force coeffcients
load(strcat('HydroCoeffs/Cylinder/Excite.mat'));  %Load excitation force co-efficients - "Fe"
load(strcat('HydroCoeffs/Cylinder/w.mat'));  %Load corresponding frequencies - "w"

%Calculate excitation force
wave_w=2*pi/Period;         %Calculate input wave frequency
ind=find(w<=wave_w);        %Find the excitation force co-efficient for that frequency
Fex_HeaveCoefficient = Fe(ind(end),3);
Fex_PitchCoefficient = Fe(ind(end),5);
Fe3=Amplitude*(real(Fex_HeaveCoefficient)*cos(wave_w*Time)+imag(Fex_HeaveCoefficient)*sin(wave_w*Time));
Fe5=Amplitude*(real(Fex_PitchCoefficient)*cos(wave_w*Time)+imag(Fex_PitchCoefficient)*sin(wave_w*Time));

Fe5=Fe5.*Ramp;
Fe3=Fe3.*Ramp;

%Excitation force vector for Simulink
Fe=[Time',Fe3',Fe5'];


%--------------------------------------------------------------------------
%%% Simulate %%%
%--------------
sim('Sim.slx')


%--------------------------------------------------------------------------
%%% Post Process %%%
%-------------------
t=simout.time;
heave=simout.Data(:,1);
pitch=simout.Data(:,2);
X3=interp1(t,heave,Time);
X5=interp1(t,pitch,Time)/6.28*360;




%Plot the input wave, heave and pitch displacement
figure, subplot(3,1,1), plot(t,wave), ylabel('Input wave (m)')
subplot(3,1,2), plot(t,X3), ylabel('Heave (m)'), grid on
subplot(3,1,3), plot(t,X5), ylabel('Pitch (deg)'), ,grid on, xlabel('Time (s)')

