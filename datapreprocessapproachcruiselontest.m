addpath('C:\Users\Chris\Documents\THESIShome\Sparse Identification\sparsedynamics\utils')
addpath('C:\Users\Chris\Documents\THESIShome\SIDPAC_ver_4.1\SIDPAC')
clear all
close all
clc
%====================================================
% Preprocess the text file so we can read it in with readtable().
fInput = fopen('Datalearapproachdoubsandsweepand3211.txt', 'rt');
fOut = fopen('temp.txt', 'wt');
str = fgetl(fInput);
% Get rid of bad characters from first line.
str(str==',') = [];
str(str=='_') = [];
lineCount = 1;
tic
while ischar(str);
  fprintf(1, 'Read line #%d\n', lineCount);
  % Get rid of bad characters.
  str(str==' ') = []; % Remove spaces
  str = strrep(str, '|', ','); % Replace bar by comma.
  % Write out good line to new file.
  fprintf(fOut, '%s\n', str);
  % Retrieve the next line.
  str = fgetl(fInput);
  lineCount = lineCount + 1;
end
toc
fclose(fInput);
fclose(fOut);

%====================================================
% Now read the table from the temporary file.
t = readtable('temp.txt', 'delimiter', ',');
% Get rid of temporary file because we don't need it anymore.
delete('temp.txt');

%  t(:,[79]) = [];
 ttest=t;
%  ttest([685:end],:) = [];
%  ttest([length(ttest.realtime)-100:length(ttest.realtime)],:) = [];

%% gather data to find times
Qt = ttest.Qrad_s;
Pt = ttest.Prad_s;
Rt = ttest.Rrad_s;
thetat = ttest.pitchdeg*pi/180;
phit = ttest.rolldeg*pi/180;
alphat = ttest.alphadeg*pi/180;
betat = ttest.betadeg*pi/180;
vtruet = ttest.Vtruektas*1.68781;
vt=betat.*vtruet;
wt=alphat.*vtruet;
Ut=sqrt(vtruet.^2-vt.^2-wt.^2);
det=ttest.elev1deg*pi/180; % de 172 +-24°
dat=(ttest.Railn1deg/2+(-1*ttest.Lailn1deg/2))*pi/180; %da 172 +-15
drt=ttest.rudr1deg*pi/180; %dr 172 +-17
dtt=ttest.thro1part;
% h=ttest.alt1ftmsl;
timet=ttest.realtime;
timetott=ttest.totltime;
timet = timet - timet(1);

% plot data
figure
subplot(5,1,1)
plot(timet,vtruet)
xlabel('time (s)')
ylabel('vtrue (ft/s)')
subplot(5,1,2)
plot(timet,alphat*180/pi)
xlabel('time (s)')
ylabel('Angle of attack (deg)')
subplot(5,1,3)
plot(timet,Qt*180/pi)
xlabel('time (s)')
ylabel('Pitch Rate (deg/s)')
subplot(5,1,4)
plot(timet,thetat*180/pi)
xlabel('time (s)')
ylabel('Pitch (deg)')
subplot(5,1,5)
plot(timet,det*180/pi)
xlabel('time (s)')
ylabel('elevator (deg)')


%% gather variables for doublet
Qdoub = ttest.Qrad_s(24200:25000,:);
Pdoub = ttest.Prad_s(24200:25000,:);
Rdoub = ttest.Rrad_s(24200:25000,:);
thetadoub = ttest.pitchdeg(24200:25000,:)*pi/180;
thetadoub=thetadoub-4*pi/180;
phidoub = ttest.rolldeg(24200:25000,:)*pi/180;
alphadoub = ttest.alphadeg(24200:25000,:)*pi/180;
alphadoub=alphadoub-4*pi/180;
betadoub = ttest.betadeg(24200:25000,:)*pi/180;
vtruedoub = ttest.Vtruektas(24200:25000,:)*1.68781;
vtruedoub = vtruedoub-170;
vdoub=betadoub.*vtruedoub;
wdoub=alphadoub.*vtruedoub;
Udoub=sqrt(vtruedoub.^2-vdoub.^2-wdoub.^2);
dedoub=ttest.elev1deg(24200:25000,:)*pi/180; % de 172 +-24°
dedoub = dedoub+2.1*pi/180;
dadoub=(ttest.Railn1deg(24200:25000,:)/2+(-1*ttest.Lailn1deg(24200:25000,:)/2))*pi/180; %da 172 +-15
drdoub=ttest.rudr1deg(24200:25000,:)*pi/180; %dr 172 +-17
dtdoub=ttest.thro1part(24200:25000,:);
% h=ttest.alt1ftmsl;
timedoub=ttest.realtime(24200:25000,:);
timedoub = timedoub-timedoub(1);

for ind = 1:length(timedoub)-1
    delt(ind) = timedoub(ind+1)-timedoub(ind);
    ind = ind+1;
end
wsample = mean(delt);
fs = round(1/wsample);

% plot data
figure
subplot(5,1,1)
plot(timedoub,vtruedoub)
xlabel('time (s)')
ylabel('vtrue (ft/s)')
subplot(5,1,2)
plot(timedoub,alphadoub*180/pi)
xlabel('time (s)')
ylabel('Angle of attack (deg)')
subplot(5,1,3)
plot(timedoub,Qdoub*180/pi)
xlabel('time (s)')
ylabel('Pitch Rate (deg/s)')
subplot(5,1,4)
plot(timedoub,thetadoub*180/pi)
xlabel('time (s)')
ylabel('Pitch (deg)')
subplot(5,1,5)
plot(timedoub,dedoub*180/pi)
xlabel('time (s)')
ylabel('elevator (deg)')

% pad signals for resampling
delvpaddoub = [ zeros(4*fs,1); vtruedoub; zeros(4*fs,1)];
alppaddoub = [ zeros(4*fs,1); alphadoub; zeros(4*fs,1)];
Qpaddoub = [ zeros(4*fs,1); Qdoub; zeros(4*fs,1)];
thetapaddoub = [ zeros(4*fs,1); thetadoub; zeros(4*fs,1)];
depaddoub = [ zeros(4*fs,1); dedoub; zeros(4*fs,1)];
statesdoub = [delvpaddoub, alppaddoub, Qpaddoub, thetapaddoub];
inputsdoub = [depaddoub];
ydoub = [statesdoub inputsdoub];
% pad the time signal
t4paddoub = 0:wsample:4;
t4paddoub = t4paddoub';
timepaddoub = [t4paddoub+timedoub(1)-t4paddoub(end)-wsample; timedoub; t4paddoub+timedoub(end)+wsample]; 
[yrsdoub,trsdoub] = resample(ydoub,timepaddoub,fs);

% plot to check resampling
figure, plot(trsdoub,yrsdoub,'x'), hold on, plot(timepaddoub,ydoub,'o')

% chop off padded data
%7621 is the last data point
yrsdoub(1:length(t4paddoub)+1,:)=[];
yrsdoub(end-length(t4paddoub)+1:end,:)=[];
trsdoub(1:length(t4paddoub)+1,:)=[];
trsdoub(end-length(t4paddoub)+1:end,:)=[];
trsdoub = trsdoub-trsdoub(1);

%% gather data for sys id

Q = ttest.Qrad_s(26847:31515);
P = ttest.Prad_s(26847:31515);
R = ttest.Rrad_s(26847:31515);
theta = ttest.pitchdeg(26847:31515)*pi/180;
theta = theta-4*pi/180;
phi = ttest.rolldeg(26847:31515)*pi/180;
alpha = ttest.alphadeg(26847:31515)*pi/180;
alpha = alpha-4*pi/180;
beta = ttest.betadeg(26847:31515)*pi/180;
vtrue = ttest.Vtruektas(26847:31515)*1.68781;
vtrue = vtrue-170;
v=beta.*vtrue;
w=alpha.*vtrue;
U=sqrt(vtrue.^2-v.^2-w.^2);
de=ttest.elev1deg(26847:31515)*pi/180; % de 172 +-24°
de = de+2.1*pi/180;
da=(ttest.Railn1deg(26847:31515)/2+(-1*ttest.Lailn1deg(26847:31515)/2))*pi/180; %da 172 +-15
dr=ttest.rudr1deg(26847:31515)*pi/180; %dr 172 +-17
dt=ttest.thro1part(26847:31515);
% h=ttest.alt1ftmsl;
time=ttest.realtime(26847:31515);
timetot=ttest.totltime(26847:31515);
time = time - time(1);


% plot data
figure
fig = gcf
fig.Units = 'inches' 
fig.Position = [3 0.5 8 7]
fig.Color = [1 1 1]
subplot(5,1,1)
plot(time,vtrue)
xlabel('time (s)')
ylabel('\Deltavtrue (ft/s)')
subplot(5,1,2)
plot(time,alpha*180/pi)
xlabel('time (s)')
ylabel('\DeltaAngle of attack (deg)')
subplot(5,1,3)
plot(time,Q*180/pi)
xlabel('time (s)')
ylabel('Pitch Rate (deg/s)')
subplot(5,1,4)
plot(time,theta*180/pi)
xlabel('time (s)')
ylabel('Pitch (deg)')
subplot(5,1,5)
plot(time,de*180/pi)
xlabel('time (s)')
ylabel('elevator (deg)')

%% SINDY
% calculate sample rate
for ind = 1:length(time)-1
    delt(ind) = time(ind+1)-time(ind);
    ind = ind+1;
end
wsample = mean(delt);
fs = round(1/wsample);

% calculate delta v
delvtrue = vtrue-mean(vtrue(1:5*fs));

% pad signals for resampling
delvpad = [ zeros(4*fs,1); vtrue; zeros(4*fs,1)];
alppad = [ zeros(4*fs,1); alpha; zeros(4*fs,1)];
Qpad = [ zeros(4*fs,1); Q; zeros(4*fs,1)];
thetapad = [ zeros(4*fs,1); theta; zeros(4*fs,1)];
depad = [ zeros(4*fs,1); de; zeros(4*fs,1)];
states = [delvpad, alppad, Qpad, thetapad];
inputs = [depad];
y = [states inputs];
% pad the time signal
t4pad = 0:wsample:4;
t4pad = t4pad';
timepad = [t4pad+time(1)-t4pad(end)-wsample; time; t4pad+time(end)+wsample]; 
[yrs,trs] = resample(y,timepad,fs);

% plot to check resampling
figure, plot(trs,yrs,'x'), hold on, plot(timepad,y,'o')

% chop off padded data
%7621 is the last data point
yrs(1:length(t4pad)+1,:)=[];
yrs(end-length(t4pad)+1:end,:)=[];
trs(1:length(t4pad)+1,:)=[];
trs(end-length(t4pad)+1:end,:)=[];
trs = trs-trs(1);
figure, plot(trs,yrs)

% find derivatives
X = yrs(:,1:4);
U = yrs(:,5);
Xdot = deriv(X,1/fs);

% 0 - longitudinal states
% 1 - lateral-directional states
% 2 - all states
modelType = 0;

% t = time1720;
% startTime = 0;
% endTime = t(end);
% startIdx = length(t) - length(t(t>startTime))+1;
% endIdx = length(t(t<endTime));

x_long = [vtrue,alpha,Q,theta];
x_latdir = [beta,P,R,phi];

u_long = [de];
u_latdir = [da,dr];

x_all = [x_long,x_latdir];
u_all = [u_long,u_latdir];

switch modelType
    case 0
        x = x_long;
        u = u_long;
        stateNames = {'Vtrue','alpha','q','theta'};
        inputNames = {'elevator'};
    case 1
        x = x_latdir;
        u = u_latdir;
        stateNames = {'beta','P','R','phi'};
        inputNames = {'aileron','rudder'};
    case 2
        x = x_all;
        u = u_all;
        stateNames = {'u','w','q','theta','v','P','R','phi'};
        inputNames = {'elevator','aileron','rudder'};
    otherwise
        x = x_all;
        u = u_all;
        stateNames = {'Vtrue','alpha','q','theta','beta','P','R','phi'};
        inputNames = {'elevator','aileron','rudder'};
end

nStates = size(x,2);
nInputs = size(u,2);

% plot
% plotstates(u,t,x)

% X = [x];
% U = [u];

Xaug = [X U];



% get derivatives of states
% h = dt; % step size
% Xdot = derivative(X, t, 2);
% % Xdot = diff(X)/h;
% 
% Xaug = Xaug(startIdx:endIdx,:);
% Xdot = Xdot(startIdx:endIdx,:);

n = nStates + nInputs;

% options for poolData
% polyorder = 1; % <=1 is linear ; >1 is nonlinear
% usesine = 0;
% lambda = 2.5e-2;
% Theta = poolData(Xaug,n,polyorder,usesine);
lambda = .1;
Theta = Xaug;
Theta_norm = zeros(size(Theta,2),1); %zeros(size(Theta,2),1);ones
for i = 1:size(Theta,2)
   Theta_norm(i) = norm(Theta(:,i));
   Theta(:,i) = Theta(:,i)./Theta_norm(i);
end

% try theta without norming
% Theta = Xaug;
% lambda = 0.1;
 Xi = sparsifyDynamics(Theta,Xdot,lambda,nStates);

for i = 1:size(Theta,2)
   Xi(i,:) = Xi(i,:)./Theta_norm(i);
end

% syms x1 x2 x3 x4 x5 x6 x7 x8 x9
% syms u1 u2 u3 u4
% states = [x1 x2 x3 x4 x5 x6 x7 x8 x9].';
% inputs = [u1 u2 u3 u4].';
% xvec = [states(1:nStates)];
% uvec = [inputs(1:nInputs)];
% 
% % create state and input variable vector
% stateinputvec = [1;states(1:nStates);inputs(1:nInputs)];
% 
% % create equations from the sparse identified system
% sysID_ode = vpa(Xi.'*stateinputvec,2);


% identified system matrices
A_sysID = Xi(1:nStates,:).';
B_sysID = Xi(nStates+1:end,:).';
C_sysID = eye(size(A_sysID));
D_sysID = zeros(size(A_sysID,1),size(B_sysID,2));

eA = eig(A_sysID);

%%

% A_sysID(1,1) = .0005454;
% A_sysID(4,1) = 5.921e-7;
% B_sysID(1,1) = 8.401;
sysID = ss(A_sysID,B_sysID,C_sysID,D_sysID,'StateName',stateNames,'InputName',inputNames)

figure('color','w')
pzmap(sysID)
box on

damp(sysID)


%% sim inputs

x0 = [X(1,1) X(1,2) X(1,3) X(1,4)];
lsimout = lsim(sysID,U,trs,x0);
figure, 
subplot(4,1,1)
plot(trs,lsimout(:,1),trs,X(:,1))
subplot(4,1,2)
plot(trs,lsimout(:,2),trs,X(:,2))
subplot(4,1,3)
plot(trs,lsimout(:,3),trs,X(:,3))
subplot(4,1,4)
plot(trs,lsimout(:,4),trs,X(:,4))


%% OEM
% p0 = [-.1;-.1;-.3;0;-1;1;0;0.001;-6;-1;0;0;0;-6;0];
% [y,p,crb,rr]=oe('oetest2',p0,U,trs,x0,0,X);
% 
% Aoe = [p(1),p(2),p(3),-32.17;
%    p(4),p(5),p(6),p(7)
%    p(8),p(9),p(10),p(11)
%    0,0,1,0];
% Boe =[p(12);p(13);p(14);p(15)];

Aoe = [-0.0713521393046585,4.10163500784009,-10.0117116823501,-32.1700000000000;-0.00232531389613736,-1.29955290911919,0.938105444125845,-0.00549924418307547;0.000967112270037410,-6.07244924095437,-1.24137205193454,-0.00176115968607703;0,0,1,0];
Boe = [-9.56187473567135;-0.198675635695217;-6.55412126564816;-0.0720622373010956];   
sysOE = ss(Aoe,Boe,C_sysID,D_sysID,'StateName',stateNames,'InputName',inputNames)

x0 = [X(1,1) X(1,2) X(1,3) X(1,4)];
lsimout = lsim(sysOE,U,trs,x0);
figure, 
subplot(4,1,1)
plot(trs,lsimout(:,1),trs,X(:,1))
subplot(4,1,2)
plot(trs,lsimout(:,2),trs,X(:,2))
subplot(4,1,3)
plot(trs,lsimout(:,3),trs,X(:,3))
subplot(4,1,4)
plot(trs,lsimout(:,4),trs,X(:,4))

%% compare
figure
pzmap(sysOE,sysID)
legend('Output Error Method','SINDY')

x0 = [yrsdoub(1,1) yrsdoub(1,2) yrsdoub(1,3) yrsdoub(1,4)];
lsimoutdoubOE = lsim(sysOE,yrsdoub(:,5),trsdoub,x0);
lsimoutdoubSINDY = lsim(sysID,yrsdoub(:,5),trsdoub,x0);
figure, 
subplot(5,1,1)
plot(trsdoub,yrsdoub(:,1),trsdoub,lsimoutdoubOE(:,1),trsdoub,lsimoutdoubSINDY(:,1))
xlabel('time (s)')
ylabel('Delta vtrue (ft/s)')
legend('X-Plane data','Output Error Method','SINDY Method')

subplot(5,1,2)
plot(trsdoub,yrsdoub(:,2)*180/pi,trsdoub,lsimoutdoubOE(:,2)*180/pi,trsdoub,lsimoutdoubSINDY(:,2)*180/pi)
xlabel('time (s)')
ylabel('Delta Angle of attack (deg)')
legend('X-Plane data','Output Error Method','SINDY Method')

subplot(5,1,3)
plot(trsdoub,yrsdoub(:,3)*180/pi,trsdoub,lsimoutdoubOE(:,3)*180/pi,trsdoub,lsimoutdoubSINDY(:,3)*180/pi)
xlabel('time (s)')
ylabel('Pitch Rate (deg/s)')
legend('X-Plane data','Output Error Method','SINDY Method')

subplot(5,1,4)
plot(trsdoub,yrsdoub(:,4)*180/pi,trsdoub,lsimoutdoubOE(:,4)*180/pi,trsdoub,lsimoutdoubSINDY(:,4)*180/pi)
xlabel('time (s)')
ylabel('Pitch Angle (deg)')
legend('X-Plane data','Output Error Method','SINDY Method')

subplot(5,1,5)
plot(trsdoub,yrsdoub(:,5)*180/pi)
xlabel('time (s)')
ylabel('Delta Elevator Input (deg)')

% Calculate
qbar = 34.3;
S = 230;
m = 13000/32.174;
u0 = 170;
CD1 = 0.256;
CL0 = 1.64;
cbar = 7;
Iyy = 18800;
rho = 0.0023;
theta0 = 0.0698;

[Coefs] = lonStabDerivs(A_sysID, B_sysID, m, u0, Iyy, rho, S, CL0, theta0, cbar)
Coefs.Properties.VariableNames = {'Sindy'}

%% Calculate coefs OE
[CoefsOE] = lonStabDerivs(Aoe, Boe, m, u0, Iyy, rho, S, CL0, theta0, cbar)
CoefsOE.Properties.VariableNames = {'Output_Error'}
%% Calculate Roskam SS Representation
qbar = 34.3;
S = 230;
m = 13000/32.174;
U1 = 170;
CD1 = 0.256;
CL1 = 1.64;
cbar = 7;
Iyy = 18800;

Ctx1 = .2560;
Cm1 = 0;
Cmt1 = 0;
Cd0 = .0431;
Cdu = 0;
Cdalp =1.06;
Ctxu = -.6;
CL0 = 1.2;
CLu = .04;
CLalp = 5.04;
CLalpdot =1.6;
CLq = 4.1;
CM0 = .047;
CMu = -.01;
CMalp = -.66;
CMalpdot =-5;
CMq = -13.5;
CMtu = 0.006;
CMTalp =0;
CDde = 0;
CLde = .4;
CMde = -.98;


udotu = qbar*S/m*(-(Cdu+2*CD1)/U1+(Ctxu+2*Ctx1)/U1);
udotalp = qbar*S/m*(-(Cdalp));
udotq = 0;
udottheta = -32.17;
udotde = -qbar*S*CDde;

alpdotu = -qbar*S/(m*U1^2)*(CLu+2*CL1);
alpdotalp = -qbar*S/(m*U1)*(CLalp+CD1)-qbar*S*cbar/(2*m*U1^2);
alpdotq = -qbar*S*cbar/(2*m*U1^2)*CLq+1;
alpdottheta = 0;
alpdotde = -qbar*S*CLde/(m*U1);

qdotu = qbar*S*cbar*(CMu+2*Cm1)/(Iyy*U1)+qbar*S*cbar*CMtu/(Iyy*U1);
qdotalp = qbar*S*cbar*CMalp/Iyy+qbar*S*cbar^2*CMalpdot/(Iyy*2*U1);
qdotq = qbar*S*cbar^2*CMq/(Iyy*2*U1);
qdottheta=0;
qdotde = qbar*S*cbar*CMde/(Iyy);

Aroskam = [udotu udotalp udotq udottheta;
            alpdotu alpdotalp alpdotq alpdottheta;
            qdotu qdotalp qdotq qdottheta;
            0 0 1 0];
Broskam = [udotde;alpdotde;qdotde;0];



[CoefsRoskam] = lonStabDerivs(Aroskam, Broskam, m, u0, Iyy, rho, S, CL1, theta0, cbar)
CoefsRoskam.Properties.VariableNames = {'Roskam'}

sysR = ss(Aroskam,Broskam,C_sysID,D_sysID,'StateName',stateNames,'InputName',inputNames);

%% Compare doublets with Roskam


x0 = [yrsdoub(1,1) yrsdoub(1,2) yrsdoub(1,3) yrsdoub(1,4)];
lsimoutdoubOE = lsim(sysOE,yrsdoub(:,5),trsdoub,x0);
lsimoutdoubSINDY = lsim(sysID,yrsdoub(:,5),trsdoub,x0);
lsimoutdoubR = lsim(sysR,yrsdoub(:,5),trsdoub,x0);
figure, 

fig = gcf
fig.Units = 'inches' 
fig.Position = [3 0.5 8 7]
fig.Color = [1 1 1]
subplot(5,1,1)
plot(trsdoub,yrsdoub(:,1),trsdoub,lsimoutdoubOE(:,1),trsdoub,lsimoutdoubSINDY(:,1),trsdoub,lsimoutdoubR(:,1))
xlabel('time (s)')
ylabel('vtrue (ft/s)')
legend('X-Plane data','Output Error Method','SINDY Method','Roskam','Position',[.85 0.89 0.1 0.1])
grid on
subplot(5,1,2)
plot(trsdoub,yrsdoub(:,2)*180/pi,trsdoub,lsimoutdoubOE(:,2)*180/pi,trsdoub,lsimoutdoubSINDY(:,2)*180/pi,trsdoub,lsimoutdoubR(:,2)*180/pi)
xlabel('time (s)')
ylabel('Angle of attack (deg)')
%legend('X-Plane data','Output Error Method','SINDY Method','Roskam')
grid on
subplot(5,1,3)
plot(trsdoub,yrsdoub(:,3)*180/pi,trsdoub,lsimoutdoubOE(:,3)*180/pi,trsdoub,lsimoutdoubSINDY(:,3)*180/pi,trsdoub,lsimoutdoubR(:,3)*180/pi)
xlabel('time (s)')
ylabel('pitch rate (deg/s)')
%legend('X-Plane data','Output Error Method','SINDY Method','Roskam')
grid on
subplot(5,1,4)
plot(trsdoub,yrsdoub(:,4)*180/pi,trsdoub,lsimoutdoubOE(:,4)*180/pi,trsdoub,lsimoutdoubSINDY(:,4)*180/pi,trsdoub,lsimoutdoubR(:,4)*180/pi)
xlabel('time (s)')
ylabel('Pitch Angle (deg)')
%legend('X-Plane data','Output Error Method','SINDY Method','Roskam')
grid on
subplot(5,1,5)
plot(trsdoub,yrsdoub(:,5)*180/pi)
xlabel('time (s)')
ylabel('Elevator Input (deg)')
grid on
%% Concatenate coef table
LonStab = [CoefsRoskam Coefs CoefsOE]
%% Compare pole - zero
figure
fig = gcf
fig.Units = 'inches' 
fig.Position = [3 2 5 5]
fig.Color = [1 1 1]
pzmap(sysID,sysOE,sysR,'k')
title('')
legend('SINDY','Output-Error','Roskam','Position',[0.7 0.75 0.1 0.1])

damp(sysID)
damp(sysOE)
damp(sysR)
%% Compare Bode Plots
% delta u
figure
fig = gcf
fig.Units = 'inches' 
fig.Position = [3 2 8 5]
fig.Color = [1 1 1]
lonbodedu = bodeplot(sysID,sysOE,sysR)

lonbodedu.OutputName = {'\Deltau','Angle of Attack','Pitch Rate','Pitch Angle'};
lonbodedu.OutputVisible = {'on','off','off','off'}
xlim([0.01 10])
title('')
grid on
legend('SINDY','Output-Error','Roskam','Position',[0.75 0.75 0.1 0.1])
lonbodedu.AxesGrid.RowLabel = ''


p = getoptions(lonbodedu); 
p.PhaseWrapping = 'on' 
setoptions(lonbodedu,p)

% delta alpha
figure
fig = gcf
fig.Units = 'inches' 
fig.Position = [3 2 8 5]
fig.Color = [1 1 1]
lonbodedu = bodeplot(sysID,sysOE,sysR)

lonbodedu.OutputName = {'\Deltau','Angle of Attack','Pitch Rate','Pitch Angle'};
lonbodedu.OutputVisible = {'off','on','off','off'}
xlim([0.01 10])
title('')
grid on
legend('SINDY','Output-Error','Roskam','Position',[0.77 0.77 0.1 0.1])
lonbodedu.AxesGrid.RowLabel = ''



% q
figure
fig = gcf
fig.Units = 'inches' 
fig.Position = [3 2 8 5]
fig.Color = [1 1 1]
lonbodedu = bodeplot(sysID,sysOE,sysR)

lonbodedu.OutputName = {'\Deltau','Angle of Attack','Pitch Rate','Pitch Angle'};
lonbodedu.OutputVisible = {'off','off','on','off'}
xlim([0.01 10])
title('')
grid on
legend('SINDY','Output-Error','Roskam','Position',[0.77 0.52 0.1 0.1])
lonbodedu.AxesGrid.RowLabel = ''


p = getoptions(lonbodedu); 
p.PhaseWrapping = 'on' 
setoptions(lonbodedu,p)

% delta theta
figure
fig = gcf
fig.Units = 'inches' 
fig.Position = [3 2 8 5]
fig.Color = [1 1 1]
lonbodedu = bodeplot(sysID,sysOE,sysR)

lonbodedu.OutputName = {'\Deltau','Angle of Attack','Pitch Rate','Pitch Angle'};
lonbodedu.OutputVisible = {'off','off','off','on'}
xlim([0.01 10])
title('')
grid on
legend('SINDY','Output-Error','Roskam','Position',[0.75 0.75 0.1 0.1])
lonbodedu.AxesGrid.RowLabel = ''


p = getoptions(lonbodedu);
p.PhaseWrapping = 'on' 
% p.PhaseMatching = 'on' 
% p.PhaseMatchingFreq = 0
% p.PhaseMatchingValue = 0
setoptions(lonbodedu,p)



















