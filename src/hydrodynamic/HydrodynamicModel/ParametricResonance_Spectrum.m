%--------------------------------------------------------------------------
%%% Constants %%%
%----------------
% Density
rho = 1003;     %[kg/m^3]
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
%Initia Positions
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

Ramp = ones(1,length(Time));
for ri=1:floor(10/dt)
    Ramp(ri)=0;
end

for ri = floor(10/dt)+1:floor(300/dt)
Ramp(ri)=Ramp(ri)-((floor(300/dt)-(ri-floor(10/dt)))/(floor(300/dt)));
end




%--------------------------------------------------------------------------
%%% Excitation force %%%
%-----------------------
%Load Excitation force coeffcients
load(strcat('HydroCoeffs/Cylinder/Excite.mat'));  %Load excitation force co-efficients - "Fe"
load(strcat('HydroCoeffs/Cylinder/w.mat'));  %Load corresponding frequencies - "w"
FE=Fe;
load('Phases.mat')
ww=w(1:300);
dw=w(2)-w(1);
S=ww;   %Initialise vector to hold Energy Spectrum
H=ww;   %Initialise vector to hold Amplitude spectrum

%--------------------------------------------------------------------------
%%% Loop through wave frequencies %%%%
percent=0.8;    %Percentage of parametric frequency
Hs=3.5;           %Significant wave height

wp=2*w5*percent; %Peak wave frequency as a percent of parametric frequency
Fex_Heave=0*Time;
Fex_Pitch=0*Time;
wave=0*Time;
 
for ii =1:1:length(ww)
   
    
    S(ii)=5/16*wp^4/ww(ii)^5*Hs^2*exp(-1.25*(wp/ww(ii))^4);   %Pierson-Moskowitz Energy Spectrum
   H(ii)=sqrt(2*S(ii)*dw);                                 %Calculate amplitude spectrum from energy spectrum
  
       Fex_HeaveCoefficient = FE(ii,3);  %Find the excitation force co-efficient for that frequency
    Fex_PitchCoefficient = FE(ii,5);
    W=ww(ii)+dw*(-0.5+w_phases(ii));

 wave = wave+H(ii)*sin(W*Time+w_phases(ii)); % Create input wave signal[m]
    
Fex_Heave=Fex_Heave+ H(ii)*(real(Fex_HeaveCoefficient)*cos(W*Time + 2*pi*w_phases(ii))+imag(Fex_HeaveCoefficient)*sin(W*Time+ 2*pi*w_phases(ii)));
Fex_Pitch=Fex_Pitch+ H(ii)*(real(Fex_PitchCoefficient)*cos(W*Time+ 2*pi*w_phases(ii))+imag(Fex_PitchCoefficient)*sin(W*Time+ 2*pi*w_phases(ii)));
end


wave=wave.*Ramp;

Wave=[Time',wave'];  %Wave vector for Simulink
%Excitation force vector for Simulink
Fex_Heave=Fex_Heave.*Ramp;
Fex_Pitch=Fex_Pitch.*Ramp;
Fe=[Time',Fex_Heave',Fex_Pitch'];






%--------------------------------------------------------------------------
%%% Simulate %%%
%--------------
sim('Sim')


%--------------------------------------------------------------------------
%%% Post Process %%%
%-------------------


t=simout.time;
pitch=simout.Data(:,2);
X3=simout.data(:,1);
X5=interp1(t,pitch,Time)/6.28*360;

%Plot the input wave, heave and pitch displacement
figure, subplot(3,1,1), plot(t,wave), ylabel('Input wave (m)')
subplot(3,1,2), plot(t,X3), ylabel('Heave (m)'), grid on
subplot(3,1,3), plot(t,X5), ylabel('Pitch (deg)'), ,grid on, xlabel('Time (s)')
