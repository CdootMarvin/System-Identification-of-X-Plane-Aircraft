% addpath('C:\Users\Chris\Documents\THESIShome\Sparse Identification\sparsedynamics\utils')
% addpath('C:\Users\Chris\Documents\THESIShome\SIDPAC_ver_4.1\SIDPAC')
addpath('C:\Users\candress\Desktop\Thesis\Data\Sparse Identification\sparsedynamics\utils')
addpath('C:\Users\candress\Desktop\Thesis\Data\SIDPAC_ver_4.1\SIDPAC')
clear all
close all
clc
%====================================================
% Preprocess the text file so we can read it in with readtable().
fInput = fopen('Dataleardoubandsweepcruise.txt', 'rt');
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

%% gather variables for doublet
Qdoub = ttest.Qrad_s(34803:35310,:);
Pdoub = ttest.Prad_s(34803:35310,:);
Rdoub = ttest.Rrad_s(34803:35310,:);
thetadoub = ttest.pitchdeg(34803:35310,:)*pi/180;
thetadoub = thetadoub-thetadoub(1);
phidoub = ttest.rolldeg(34803:35310,:)*pi/180;
alphadoub = ttest.alphadeg(34803:35310,:)*pi/180;
alphadoub = alphadoub-alphadoub(1);
betadoub = ttest.betadeg(34803:35310,:)*pi/180;
vtruedoub = ttest.Vtruektas(34803:35310,:)*1.68781;
vtruedoub = vtruedoub-vtruedoub(55);
vdoub=betadoub.*vtruedoub;
wdoub=alphadoub.*vtruedoub;
Udoub=sqrt(vtruedoub.^2-vdoub.^2-wdoub.^2);
dedoub=ttest.elev1deg(34803:35310,:)*pi/180; % de 172 +-24°
dedoub = dedoub-dedoub(1);
dadoub=(ttest.Railn1deg(34803:35310,:)/2+(-1*ttest.Lailn1deg(34803:35310,:)/2))*pi/180; %da 172 +-15
drdoub=ttest.rudr1deg(34803:35310,:)*pi/180; %dr 172 +-17
dtdoub=ttest.thro1part(34803:35310,:);
% h=ttest.alt1ftmsl;
timedoub=ttest.realtime(34803:35310,:);
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
t4paddoub = 0:wsample:4-wsample;
t4paddoub = t4paddoub';
timepaddoub = [t4paddoub+timedoub(1)-t4paddoub(end)-wsample; timedoub; t4paddoub+timedoub(end)+wsample]; 
[yrsdoub,trsdoub] = resample(ydoub,timepaddoub,fs);

% plot to check resampling


% chop off padded data
%7621 is the last data point
yrsdoub(1:length(t4paddoub)+1,:)=[];
yrsdoub(end-length(t4paddoub)+1:end,:)=[];
trsdoub(1:length(t4paddoub)+1,:)=[];
trsdoub(end-length(t4paddoub)+1:end,:)=[];
trsdoub = trsdoub-trsdoub(1);

figure, plot(trsdoub,yrsdoub,'x')

% gather data for sys id

Q = ttest.Qrad_s(39305:45243);
P = ttest.Prad_s(39305:45243);
R = ttest.Rrad_s(39305:45243);
theta = ttest.pitchdeg(39305:45243)*pi/180;
phi = ttest.rolldeg(39305:45243)*pi/180;
alpha = ttest.alphadeg(39305:45243)*pi/180-ttest.alphadeg(34803)*pi/180;
beta = ttest.betadeg(39305:45243)*pi/180;
vtrue = ttest.Vtruektas(39305:45243)*1.68781-662.446;
v=beta.*vtrue;
w=alpha.*vtrue;
U=sqrt(vtrue.^2-v.^2-w.^2);
de=ttest.elev1deg(39305:45243)*pi/180;% de 172 +-24°
de = de-de(1);
da=(ttest.Railn1deg(39305:45243)/2+(-1*ttest.Lailn1deg(39305:45243)/2))*pi/180; %da 172 +-15
dr=ttest.rudr1deg(39305:45243)*pi/180; %dr 172 +-17
dt=ttest.thro1part(39305:45243);
% h=ttest.alt1ftmsl;
time=ttest.realtime(39305:45243);
timetot=ttest.totltime(39305:45243);
time = time - time(1);


% plot data
figure
fig = gcf
fig.Units = 'inches' 
fig.Position = [3 0.5 8 7]
fig.Color = [1 1 1]
subplot(5,1,1)
plot(time,vtrue)
xlabel('Time (s)')
ylabel('\Deltavtrue (ft/s)')
subplot(5,1,2)
plot(time,alpha*180/pi)
xlabel('Time (s)')
ylabel('\DeltaAngle of attack (deg)')
subplot(5,1,3)
plot(time,Q*180/pi)
xlabel('Time (s)')
ylabel('Pitch Rate (deg/s)')
subplot(5,1,4)
plot(time,theta*180/pi)
xlabel('Time (s)')
ylabel('Pitch Angle (deg)')
subplot(5,1,5)
plot(time,de*180/pi)
xlabel('Time (s)')
ylabel('Elevator (deg)')


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
% p0 = [0;-15;-1.44;0;-1;1.1;0;0.001;-21.6;-1;0;0;0;-24.7;0];
% [y,p,crb,rr]=oe('oetest2',p0,U,trs,x0,0,X);
% 
% Aoe = [p(1),p(2),p(3),-32.17;
%    p(4),p(5),p(6),p(7)
%    p(8),p(9),p(10),p(11)
%    0,0,1,0];
% Boe =[p(12);p(13);p(14);p(15)];

Aoe = [-0.0252867199403407,-15.0268519034650,-1.47565924675249,-32.1700000000000;-0.000113872258360225,-1.92057012848856,0.998403891062329,0.0104854206377752;0.000969870093508476,-26.3462715672886,-1.13257568109726,-0.0109256107316613;0,0,1,0];
Boe = [0.0248686977262801;-0.382662084175575;-30.9431409412143;-0.0369598938018110];

sysOE = ss(Aoe,Boe,C_sysID,D_sysID,'StateName',stateNames,'InputName',inputNames);

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

x0 = [yrsdoub(1,1) yrsdoub(1,2) yrsdoub(1,3) yrsdoub(1,4)];
lsimoutdoubOE = lsim(sysOE,yrsdoub(:,5),trsdoub,x0);
lsimoutdoubSINDY = lsim(sysID,yrsdoub(:,5),trsdoub,x0);
figure, 
subplot(5,1,1)
plot(trsdoub,yrsdoub(:,1),trsdoub,lsimoutdoubOE(:,1),trsdoub,lsimoutdoubSINDY(:,1))
xlabel('time (s)')
ylabel('\Deltavtrue (ft/s)')
legend('X-Plane data','Output Error Method','SINDY Method')
subplot(5,1,2)
plot(trsdoub,yrsdoub(:,2)*180/pi,trsdoub,lsimoutdoubOE(:,2)*180/pi,trsdoub,lsimoutdoubSINDY(:,2)*180/pi)
xlabel('time (s)')
ylabel('\DeltaAngle of attack (deg)')
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
ylabel('Elevator Input (deg)')

% Calculate
qbar = 134.6;
S = 230;
m = 13000/32.174;
u0 = 665;
CD1 = 0.0279;
CL0 = 0.28;
cbar = 7;
Iyy = 18800;
rho = qbar*2/677^2;
theta0 = 0;

[Coefs] = lonStabDerivs(A_sysID, B_sysID, m, u0, Iyy, rho, S, CL0, theta0, cbar)
Coefs.Properties.VariableNames = {'Sindy'}

%% Calculate coefs OE
[CoefsOE] = lonStabDerivs(Aoe, Boe, m, u0, Iyy, rho, S, CL0, theta0, cbar)
CoefsOE.Properties.VariableNames = {'Output_Error'}
%% Calculate Roskam SS Representation
qbar = 134.6;
S = 230;
m = 13000/32.174;
U1 = 677;
CD1 = 0.0279;
CL1 = 0.28;
cbar = 7;
Iyy = 18800;

Ctx1 = .0279;
Cm1 = 0;
Cmt1 = 0;
Cd0 = .0216;
Cdu = 0.104;
Cdalp =.3;
Ctxu = -.07;
CL0 = .13;
CLu = .28;
CLalp = 5.84;
CLalpdot =2.2;
CLq = 4.7;
CM0 = .05;
CMu = .07;
CMalp = -.64;
CMalpdot =-6.7;
CMq = -15.5;
CMtu = -.003;
CMTalp =0;
CDde = 0;
CLde = .46;
CMde = -1.24;


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
writetable(LonStab,'LearCruiseLonStabDerivs.xlsx','WriteRowNames',true)
%% Compare pole - zero
figure
fig = gcf
fig.Units = 'inches' 
fig.Position = [3 2 8 5]
fig.Color = [1 1 1]
pzmap(sysID,sysOE,sysR,'k')
title('')
legend('SINDY','Output-Error','Roskam','Position',[0.75 0.75 0.1 0.1])

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
legend('SINDY','Output-Error','Roskam','Position',[0.6 0.52 0.1 0.1])
lonbodedu.AxesGrid.RowLabel = ''

lonbodedu.AxesGrid.YLabel = {'Phase','Magnitude'}
lonbodedu.AxesGrid.YUnits = {'Deg','DB'}

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

%%
sys2 = sysID;
sys2.A = [-0.1   0     0   -32;
           0       -1.6    1    0;
           0       -6  -6   -.25;
           0         0   0.9981    0];
damp(sys2)
figure, pzmap(sys2)
figure, step(sys2,10)
           

%% eig assign

s = tf('s');
wn = .05;
dampr = .707;
tfph = wn^2/(s^2+2*dampr*wn*s+wn^2)
damp(tfph)
DesEigVal=[-5+1i*3 -5-1i*3 -.05+1i*.05 -.05-1i*.05];
       
DesEigVec=[1,           nan,        0,          0;
           0,           0,          nan,        1;
           0,           0,          1,          nan;
           nan,         nan,        0,          0
           ];
       
ConsMat = nan(1,4);

Gain = EigStructDesign(sysID.A, sysID.B, sysID.C, DesEigVal, DesEigVec, ConsMat);

sysclose = feedback(sysID,Gain)

[veccl, eigcl] = eig(sysclose.A)


figure, pzmap(sysclose)
figure, step(sysID, sysclose,100)

%% imf
% desired short period mode
Asp = [-2   1;
       -13  -6];
Bsp = sysID.B(2:3,1);
Csp = eye(2);
Dsp = zeros(2,1);

syssp = ss(Asp,Bsp,Csp,Dsp)
damp(syssp)
syssp.OutputName = {'alpha','pitch rate'};
syssp.StateName = {'alpha','pitch rate'};
syssp.InputName = {'Elevator'};

figure, step(syssp)
figure, pzmap(syssp(1,1),syssp(2,1))

% desired phugoid mode
% assume speed and pitch rate equations, aka constant alpha
wnph = .05;
dampph = .707;
tfph = wnph^2/(s^2+2*dampph*wnph*s+wnph^2)
Aph = [-.0707   -32.2;
       .000077  0];
% A11 = Xu and A21 = Zu
Bph = sysID.B([1 4],1);
Cph = eye(2);
Dph = zeros(2,1);

sysph = ss(Aph,Bph,Cph,Dph)
damp(sysph)
syssp.OutputName = {'speed','pitch angle'};
syssp.StateName = {'speed','pitch angle'};
syssp.InputName = {'Elevator'};

figure, step(sysph)
figure, pzmap(sysph(1,1),sysph(2,1))
impulse(sysph(1,1),sysph(2,1),sysID(1,1),sysID(4,1),30)

% this determines Zu = -.000077 and Xu = -0.0707


% combine desired A matrices
Ades = [-0.0707 0 0 -32.2;
        -0.000077 -2 1 0;
        0         -13 -6 0;
        0          0  1 0];
sysdes = ss(Ades,sysID.B, sysID.C,sysID.D)   
    
figure, step(sysdes)
       
Am = Ades;
A = sysID.A;
Bprime = sysID.B;     
Qprime=eye(4);
Qprime(1,1)=1/1^2;
Qprime(2,2)=1/(1*pi/180)^2;
Qprime(3,3)=1/(1*pi/180)^2;
Qprime(4,4)=1/(1*pi/180)^2;


Q=(A-Am)'*Qprime*(A-Am);
R=(Bprime'*Qprime*Bprime+1.*eye(1));
N=(A-Am)'*Qprime*Bprime;
K=lqr(A,Bprime,Q,R,N);
syscloseimf = feedback(sysID,K)
damp(syscloseimf)
[vecclimf, eigclimf] = eig(syscloseimf.A)
[vecol, eigol] = eig(sysID.A)
figure, pzmap(syscloseimf,sysID)
figure, step(sysID, syscloseimf, sysclose, 100)
figure, step(sysID, syscloseimf,sysdes,10)

figure, impulse(sysID,syscloseimf,sysdes,10)
figure, pzmap(syscloseimf,sysID)

sysIDtf = tf(sysID);

figure, rlocus(sysIDtf(4))

Kcheck = K
Kcheck(1) = 0
fdbkcheck = feedback(sysID,Kcheck)
fdbkcheck.A-syscloseimf.A
figure, pzmap(syscloseimf,fdbkcheck)
figure, impulse(syscloseimf,fdbkcheck,100)
%% calculate parameters
t = 0:.01:6;
figure, step(syscloseimf(4),syscloseimf(4) - syscloseimf(2),t)
% gamma(3.38) = -5.75; theta(2.7) = -5.75;
Ttheta2cl = 3.38-2.7;
V = 670.5;
g = 32.174;
nalpcl = V/g*1/Ttheta2cl;
wnspcl = 4.72;
dampspcl = 0.806;
CAPcl = wnspcl^2/nalpcl

%% flight path attitude lag

t = 0:.01:6;
figure, step(sysID(4),sysID(4) - sysID(2),t)
% gamma(3.41) = -5.52; theta(2.75) = -5.52;
Ttheta2 = 3.41-2.75;
V = 670.5;
g = 32.174;
nqalp = V/g*1/Ttheta2;
wnsp = 5.17;
dampsp = 0.269; 
CAP = wnsp^2/nqalp

%% wnsp vs n/alp
% vertices
% top level 3
xwnvnalp = [1 1 100 100];
ytlvl3 = [3 100 100 30];
ytlvl2 = [2 3 30 20];
ylvl1 = [0.3 2 20 3];
yllvl2 = [.2 .3 3 2];
yllvl3 = [.1 .2 2 .1];
figure,
loglog(1,1)
hold on
fill(xwnvnalp,ytlvl3,'r')
fill(xwnvnalp,ytlvl2,'m')
fill(xwnvnalp,ylvl1,'b')
fill(xwnvnalp,yllvl2,'m')
fill(xwnvnalp,yllvl3,'r')
xlim([1 100])
ylim([0.1 100])
loglog(nalpcl,wnspcl,'ys','MarkerSize',10,'MarkerFaceColor','y','MarkerEdgeColor','k','LineWidth',2)
loglog(nqalp,wnsp,'g^','MarkerSize',10,'MarkerFaceColor','g','MarkerEdgeColor','k','LineWidth',2)
ax = gca
ax.YGrid = 'on';
ax.XGrid = 'on';
ax.Layer = 'top';
ax.GridAlpha = 0.5;
xlabel('n/alpha')
ylabel('wnsp')

%% CAP vs damping ratio
xcapvdr3 = [0.1 0.1 10 10];
xcapvdr2 = [0.2 0.2 2 2];
xcapvdr1 = [0.3 0.3 2 2];
ylvl3 = [0.01 10 10 0.01];
ylvl2 = [.038 10 10 .038];
ylvl1 = [0.085 3.6 3.6 0.085];
figure,
loglog(1,1)
hold on
fill(xcapvdr3,ylvl3,'r')
fill(xcapvdr2,ylvl2,'m')
fill(xcapvdr1,ylvl1,'b')
loglog(dampspcl,CAPcl,'ys','MarkerSize',10,'MarkerFaceColor','y','MarkerEdgeColor','k','LineWidth',2)
loglog(dampsp,CAPcl,'g^','MarkerSize',10,'MarkerFaceColor','g','MarkerEdgeColor','k','LineWidth',2)
xlim([.1 10])
ylim([0.01 10])
ax = gca
ax.YGrid = 'on';
ax.XGrid = 'on';
ax.Layer = 'top';
ax.GridAlpha = 0.5;
xlabel('Short Period Damping Ratio')
ylabel('Control Anticipation Parameter')

%% dropback
qpkpqsscl = 5.9/3.2;
DBthpqsscl = 0.8/3.2;

qpkqssba = 9/2.7;
DBthpqssba = 1.5/2.7;
xdropback = [0 0 2 2];
ylvl1 = [0 4 4 0];
ylvl2 = [3 4 4 1.625];
ylvl3 = [3.5 4 4 2.125];
figure,
plot(1,1)
hold on
fill(xdropback,ylvl1,'b')
fill(xdropback,ylvl2,'m')
fill(xdropback,ylvl3,'r')
xlim([0 2])
ylim([1 4])
ax = gca
ax.YGrid = 'on';
ax.XGrid = 'on';
ax.Layer = 'top';
ax.GridAlpha = 1;
ax.GridLineStyle = ':';
xticks(0:0.5:2)
yticks(1:0.5:4)
plot(DBthpqsscl,qpkpqsscl,'ys','MarkerSize',10,'MarkerFaceColor','y','MarkerEdgeColor','k','LineWidth',2)
loglog(DBthpqssba,qpkqssba,'g^','MarkerSize',10,'MarkerFaceColor','g','MarkerEdgeColor','k','LineWidth',2)
xlabel('Theta Dropback per Steady State Pitch Rate')
ylabel('Peak Pitch rate per Steady State Pitch Rate')


%% Pitch att response
t = 0:0.02:10;
inputcl = ones(1,length(t));
inputcl(1:5/.02) = -0.03;
inputcl(5/.02+1:end) = -.01;
cltheta = lsim(syscloseimf(4,1),inputcl,t);

inputba = ones(1,length(t));
inputba(1:5/.02) = -0.0315;
inputba(5/.02+1:end) = -.008;

batheta = lsim(sysID(4,1),inputba,t);
figure
plot(t,cltheta,'LineWidth',2), hold on
plot(t,batheta,'LineWidth',2)
grid on

xlabel('time (s)')
ylabel('Pitch Attitude Response (deg)')
legend('Closed Loop Response','Bare Airframe')

