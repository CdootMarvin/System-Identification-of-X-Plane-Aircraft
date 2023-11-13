addpath('C:\Users\Chris\Documents\THESIShome\Sparse Identification\sparsedynamics\utils')
addpath('C:\Users\Chris\Documents\THESIShome\SIDPAC_ver_4.1\SIDPAC')
clear all 
close all
clc

%====================================================
% Preprocess the text file so we can read it in with readtable().
fInput = fopen('Data172100ktslon.txt', 'rt');
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
Qdoub = ttest.Qrad_s(53075:53700,:);
Pdoub = ttest.Prad_s(53075:53700,:);
Rdoub = ttest.Rrad_s(53075:53700,:);
thetadoub = ttest.pitchdeg(53075:53700,:)*pi/180;
thetadoub = thetadoub-.44*pi/180;
phidoub = ttest.rolldeg(53075:53700,:)*pi/180;
alphadoub = ttest.alphadeg(53075:53700,:)*pi/180;
alphadoub = alphadoub-.2*pi/180;
betadoub = ttest.betadeg(53075:53700,:)*pi/180;
vtruedoub = ttest.Vtruektas(53075:53700,:)*1.68781;
vtruedoub = vtruedoub-172;
vdoub=betadoub.*vtruedoub;
wdoub=alphadoub.*vtruedoub;
Udoub=sqrt(vtruedoub.^2-vdoub.^2-wdoub.^2);
dedoub=ttest.elev1deg(53075:53700,:)*pi/180; % de 172 +-24°
dedoub=dedoub-1.029*pi/180;
dadoub=(ttest.Railn1deg(53075:53700,:)/2+(-1*ttest.Lailn1deg(53075:53700,:)/2))*pi/180; %da 172 +-15
drdoub=ttest.rudr1deg(53075:53700,:)*pi/180; %dr 172 +-17
dtdoub=ttest.thro1part(53075:53700,:);
% h=ttest.alt1ftmsl;
timedoub=ttest.realtime(53075:53700,:);
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
yrsdoub(1:length(t4paddoub)+2,:)=[];
yrsdoub(end-length(t4paddoub)-1:end,:)=[];
trsdoub(1:length(t4paddoub)+2,:)=[];
trsdoub(end-length(t4paddoub)-1:end,:)=[];
trsdoub = trsdoub-trsdoub(1);

% gather data for sys id

Q = ttest.Qrad_s(39100:42960);
P = ttest.Prad_s(39100:42960);
R = ttest.Rrad_s(39100:42960);
theta = ttest.pitchdeg(39100:42960)*pi/180-.44*pi/180;
phi = ttest.rolldeg(39100:42960)*pi/180;
alpha = ttest.alphadeg(39100:42960)*pi/180-.2*pi/180;
beta = ttest.betadeg(39100:42960)*pi/180;
vtrue = ttest.Vtruektas(39100:42960)*1.68781-172;
v=beta.*vtrue;
w=alpha.*vtrue;
U=sqrt(vtrue.^2-v.^2-w.^2);
de=ttest.elev1deg(39100:42960)*pi/180; % de 172 +-24°
de = de-1.029*pi/180;
da=(ttest.Railn1deg(39100:42960)/2+(-1*ttest.Lailn1deg(39100:42960)/2))*pi/180; %da 172 +-15
dr=ttest.rudr1deg(39100:42960)*pi/180; %dr 172 +-17
dt=ttest.thro1part(39100:42960);
% h=ttest.alt1ftmsl;
time=ttest.realtime(39100:42960);
timetot=ttest.totltime(39100:42960);
time = time - time(1);


% plot data
figure
subplot(5,1,1)
plot(time,vtrue)
xlabel('Time (s)')
ylabel('\Delta vtrue (ft/s)')
xlim([0 61])
subplot(5,1,2)
plot(time,alpha*180/pi)
xlabel('Time (s)')
ylabel({'\Delta Angle';'of attack (deg)'})
xlim([0 61])
subplot(5,1,3)
plot(time,Q*180/pi)
xlabel('Time (s)')
ylabel({'Pitch Rate';'(deg/s)'})
xlim([0 61])
subplot(5,1,4)
plot(time,theta*180/pi)
xlabel('Time (s)')
ylabel('Pitch (deg)')
xlim([0 61])
subplot(5,1,5)
plot(time,de*180/pi)
xlabel('Time (s)')
ylabel('Elevator (deg)')
xlim([0 61])


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
yrs(end-length(t4pad)-1:end,:)=[];
trs(1:length(t4pad)+1,:)=[];
trs(end-length(t4pad)-1:end,:)=[];
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
% p0 = [-.1;8;0;-32;-0.0001;-2;1;0;0.001;-20;-2;0;0;0;1;0;-.4;-.2;-20;0];
% p0 = [-.1;8;0;-0.0001;-2;1;0;0.001;-20;-2;0;-.4;-.2;-20;0];
% [y,p,crb,rr]=oe('oetest2',p0,U,trs,x0,0,X);
% 
% Aoe = [p(1),p(2),p(3),-32.17;
%    p(4),p(5),p(6),p(7)
%    p(8),p(9),p(10),p(11)
%    0,0,1,0];
% Boe =[p(12);p(13);p(14);p(15)];
Aoe = [0.0131732737110857,30.2357830045822,-32.5728008283775,-32.1700000000000;-0.00204932323803212,-3.05325039419798,0.922435175977956,-0.00246933464400217;-0.000102585055685006,-19.1853792563133,-2.58855121148555,0.0144673505464689;0,0,1,0];
Boe = [-69.5574378874566;-0.553458601924820;-26.1005268870821;0.0500161189343861];

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
plot(trsdoub,yrsdoub(:,1),'LineWidth',2),hold on, plot(trsdoub,lsimoutdoubOE(:,1),'LineWidth',2),plot(trsdoub,lsimoutdoubSINDY(:,1))
xlabel('Time (s)')
ylabel('\Delta vtrue (ft/s)')
legend('X-Plane data','Output Error Method','SINDY Method')
xlim([0 9.5])

subplot(5,1,2)
plot(trsdoub,yrsdoub(:,2)*180/pi,'LineWidth',2), hold on, plot(trsdoub,lsimoutdoubOE(:,2)*180/pi,'LineWidth',2),hold on, plot(trsdoub,lsimoutdoubSINDY(:,2)*180/pi,'LineWidth',2)
xlabel('Time (s)')
ylabel({'\Delta Angle';'of attack (deg)'})
xlim([0 9.5])

subplot(5,1,3)
plot(trsdoub,yrsdoub(:,3)*180/pi,'LineWidth',2),hold on, plot(trsdoub,lsimoutdoubOE(:,3)*180/pi,'LineWidth',2),plot(trsdoub,lsimoutdoubSINDY(:,3)*180/pi,'LineWidth',2)
xlabel('Time (s)')
ylabel({'Pitch Rate';'(deg/s)'})
xlim([0 9.5])

subplot(5,1,4)
plot(trsdoub,yrsdoub(:,4)*180/pi,'LineWidth',2),hold on, plot(trsdoub,lsimoutdoubOE(:,4)*180/pi,'LineWidth',2),plot(trsdoub,lsimoutdoubSINDY(:,4)*180/pi,'LineWidth',2)
xlabel('Time (s)')
ylabel('Pitch Angle (deg)')
xlim([0 9.5])

subplot(5,1,5)
plot(trsdoub,yrsdoub(:,5)*180/pi,'LineWidth',2)
xlabel('Time (s)')
ylabel('Elevator (deg)')
xlim([0 9.5])

%% NASA
rho = 0.0023;
U1 = 170;
u0=U1;
qbar = 1/2*rho*U1^2;
S = 174;
m = 1850/32.174;
theta0=0;

Ixx = 1030;
Iyy = 1090;
Izz = 1890;
Ixz = 90;

b = 36;
cbar = 4.9;

CD1 = 0.0;
CL1 = 0.36;


Ctx1 = .0;
Cm1 = 0;
Cmt1 = 0;
Cd0 = 0.0;
Cdu = 0;
Cdalp =.54;
Ctxu = 0;
CL0 = .36;
CLu = 0;
CLalp = 5.225;
CLalpdot =0;
CLq = 8.565;
CM0 = 0;
CMu = 0;
CMalp = -.678;
CMalpdot =-4;
CMq = -15.345;
CMtu = 0;
CMTalp =0;
CDde = 0;
CLde = .402;
CMde = -1.18;


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

%% Calculate Derivatives Sindy
[Coefs] = lonStabDerivs(sysID.A, sysID.B, m, u0, Iyy, rho, S, CL1, theta0, cbar)
Coefs.Properties.VariableNames = {'Sindy'}
%% Derivatives OE
[CoefsOE] = lonStabDerivs(sysOE.A, sysOE.B, m, u0, Iyy, rho, S, CL1, theta0, cbar)
CoefsOE.Properties.VariableNames = {'Output_Error'}

%% Concatenate coef table
LonStab = [CoefsRoskam Coefs CoefsOE]
writetable(LonStab,'SkyHawk100KtsLonStabDerivs.xlsx','WriteRowNames',true)

%% plots
% Pzmap
figure
pzmap(sysR,sysID,sysOE)
legend('NASA','Sindy','Output Error')
grid on
title('')

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
legend('SINDY','Output-Error','NASA','Position',[0.75 0.75 0.1 0.1])
lonbodedu.AxesGrid.RowLabel = ''


p = getoptions(lonbodedu); 
p.PhaseWrapping = 'on' 
setoptions(lonbodedu,p)
lonbodedu.AxesGrid.YLabel = {'Phase','Magnitude'}
lonbodedu.AxesGrid.YUnits = {'Deg','DB'}

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
legend('SINDY','Output-Error','NASA','Position',[0.6 0.52 0.1 0.1])
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
legend('SINDY','Output-Error','NASA','Position',[0.77 0.52 0.1 0.1])
lonbodedu.AxesGrid.RowLabel = ''


p = getoptions(lonbodedu); 
p.PhaseWrapping = 'on' 
setoptions(lonbodedu,p) 

lonbodedu.AxesGrid.YLabel = {'Phase','Magnitude'}
lonbodedu.AxesGrid.YUnits = {'Deg','DB'}



















