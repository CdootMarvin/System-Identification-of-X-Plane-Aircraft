addpath('C:\Users\candress\Desktop\Thesis\Data\Sparse Identification\sparsedynamics\utils')
addpath('C:\Users\candress\Desktop\Thesis\Data\SIDPAC_ver_4.1\SIDPAC')
clear all
close all
clc
%====================================================
% Preprocess the text file so we can read it in with readtable().
fInput = fopen('DataLearlatdirapproachdoubandsweep.txt', 'rt');
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

%% gather variables to find times
% gather data for sys id

Q = ttest.Qrad_s;
P = ttest.Prad_s;
R = ttest.Rrad_s;
theta = ttest.pitchdeg*pi/180-.44*pi/180;
phi = ttest.rolldeg*pi/180;
alpha = ttest.alphadeg*pi/180-.2*pi/180;
beta = ttest.betadeg*pi/180;
vtrue = ttest.Vtruektas*1.68781-172;
v=beta.*vtrue;
w=alpha.*vtrue;
U=sqrt(vtrue.^2-v.^2-w.^2);
de=ttest.elev1deg*pi/180; % de 172 +-24�
de = de-1.029*pi/180;
da=(ttest.Railn1deg/2+(-1*ttest.Lailn1deg/2))*pi/180; %da 172 +-15
dr=ttest.rudr1deg*pi/180; %dr 172 +-17
dt=ttest.thro1part;
% h=ttest.alt1ftmsl;
time=ttest.realtime;
timetot=ttest.totltime;
time = time - time(1);


% plot data
figure
subplot(5,1,1)
plot(time,beta)
xlabel('time (s)')
ylabel('beta (ft/s)')
subplot(5,1,2)
plot(time,P*180/pi)
xlabel('time (s)')
ylabel('Roll Rate (deg/s)')
subplot(5,1,3)
plot(time,R*180/pi)
xlabel('time (s)')
ylabel('Yaw Rate (deg/s)')
subplot(5,1,4)
plot(time,phi*180/pi)
xlabel('time (s)')
ylabel('Roll Angle (deg)')
subplot(5,1,5)
plot(time,da*180/pi,time,dr*180/pi)
xlabel('time (s)')
ylabel('Surface Deflection (deg)')
legend('Aileron','Rudder')

%% gather variables for aileron doublet
Pdoubail = ttest.Prad_s(17800:18450,:);
Rdoubail = ttest.Rrad_s(17800:18450,:);
phidoubail = ttest.rolldeg(17800:18450,:)*pi/180;
betadoubail = ttest.betadeg(17800:18450,:)*pi/180;
dadoubail=(ttest.Railn1deg(17800:18450,:)/2+(-1*ttest.Lailn1deg(17800:18450,:)/2))*pi/180; %da 172 +-15
drdoubail=ttest.rudr1deg(17800:18450,:)*pi/180; %dr 172 +-17
dtdoubail=ttest.thro1part(17800:18450,:);
% h=ttest.alt1ftmsl;
timedoubail=ttest.realtime(17800:18450,:);
timedoubail = timedoubail-timedoubail(1);

for ind = 1:length(timedoubail)-1
    delt(ind) = timedoubail(ind+1)-timedoubail(ind);
    ind = ind+1;
end
wsample = mean(delt);
fs = round(1/wsample);

% plot data
figure
subplot(5,1,1)
plot(timedoubail,betadoubail*180/pi)
xlabel('time (s)')
ylabel('beta (deg)')
subplot(5,1,2)
plot(timedoubail,Pdoubail*180/pi)
xlabel('time (s)')
ylabel('Roll Rate (deg/s)')
subplot(5,1,3)
plot(timedoubail,Rdoubail*180/pi)
xlabel('time (s)')
ylabel('Yaw Rate (deg/s)')
subplot(5,1,4)
plot(timedoubail,phidoubail*180/pi)
xlabel('time (s)')
ylabel('Roll Angle (deg)')
subplot(5,1,5)
plot(timedoubail,dadoubail*180/pi,timedoubail,drdoubail*180/pi)
xlabel('time (s)')
ylabel('Surface Deflection (deg)')
legend('Aileron','Rudder')

% pad signals for resampling
delvpaddoubail = [ zeros(4*fs,1); betadoubail; zeros(4*fs,1)];
alppaddoubail = [ zeros(4*fs,1); Pdoubail; zeros(4*fs,1)];
Qpaddoubail = [ zeros(4*fs,1); Rdoubail; zeros(4*fs,1)];
thetapaddoubail = [ zeros(4*fs,1); phidoubail; zeros(4*fs,1)];
depaddoubail = [ zeros(4*fs,1); dadoubail; zeros(4*fs,1)];
drpaddoubail = [ zeros(4*fs,1); drdoubail; zeros(4*fs,1)];
statesdoubail = [delvpaddoubail, alppaddoubail, Qpaddoubail, thetapaddoubail];
inputsdoubail = [depaddoubail, drpaddoubail];
ydoubail = [statesdoubail inputsdoubail];
% pad the time signal
t4paddoubail = 0:wsample:4;
t4paddoubail = t4paddoubail';
timepaddoubail = [t4paddoubail+timedoubail(1)-t4paddoubail(end)-wsample; timedoubail; t4paddoubail+timedoubail(end)+wsample]; 
[yrsdoubail,trsdoubail] = resample(ydoubail,timepaddoubail,fs);

% plot to check resampling
figure, plot(trsdoubail,yrsdoubail,'x'), hold on, plot(timepaddoubail,ydoubail,'o')

% chop off padded data
%7621 is the last data point
yrsdoubail(1:length(t4paddoubail)+2,:)=[];
yrsdoubail(end-length(t4paddoubail)-1:end,:)=[];
trsdoubail(1:length(t4paddoubail)+2,:)=[];
trsdoubail(end-length(t4paddoubail)-1:end,:)=[];
trsdoubail = trsdoubail-trsdoubail(1);
%% Gather data for rudder doublet

Pdoubrud = ttest.Prad_s(22550:23300,:);
Rdoubrud = ttest.Rrad_s(22550:23300,:);
phidoubrud = ttest.rolldeg(22550:23300,:)*pi/180;
betadoubrud = ttest.betadeg(22550:23300,:)*pi/180;
dadoubrud=(ttest.Railn1deg(22550:23300,:)/2+(-1*ttest.Lailn1deg(22550:23300,:)/2))*pi/180; %da 172 +-15
drdoubrud=ttest.rudr1deg(22550:23300,:)*pi/180; %dr 172 +-17
dtdoubrud=ttest.thro1part(22550:23300,:);
% h=ttest.alt1ftmsl;
timedoubrud=ttest.realtime(22550:23300,:);
timedoubrud = timedoubrud-timedoubrud(1);

for ind = 1:length(timedoubrud)-1
    delt(ind) = timedoubrud(ind+1)-timedoubrud(ind);
    ind = ind+1;
end
wsample = mean(delt);
fs = round(1/wsample);

% plot data
figure
subplot(5,1,1)
plot(timedoubrud,betadoubrud*180/pi)
xlabel('time (s)')
ylabel('beta (deg)')
subplot(5,1,2)
plot(timedoubrud,Pdoubrud*180/pi)
xlabel('time (s)')
ylabel('Roll Rate (deg/s)')
subplot(5,1,3)
plot(timedoubrud,Rdoubrud*180/pi)
xlabel('time (s)')
ylabel('Yaw Rate (deg/s)')
subplot(5,1,4)
plot(timedoubrud,phidoubrud*180/pi)
xlabel('time (s)')
ylabel('Roll Angle (deg)')
subplot(5,1,5)
plot(timedoubrud,dadoubrud*180/pi,timedoubrud,drdoubrud*180/pi)
xlabel('time (s)')
ylabel('Surface Deflection (deg)')
legend('Aileron','Rudder')

% pad signals for resampling
delvpaddoubrud = [ zeros(4*fs,1); betadoubrud; zeros(4*fs,1)];
alppaddoubrud = [ zeros(4*fs,1); Pdoubrud; zeros(4*fs,1)];
Qpaddoubrud = [ zeros(4*fs,1); Rdoubrud; zeros(4*fs,1)];
thetapaddoubrud = [ zeros(4*fs,1); phidoubrud; zeros(4*fs,1)];
depaddoubrud = [ zeros(4*fs,1); dadoubrud; zeros(4*fs,1)];
drpaddoubrud = [ zeros(4*fs,1); drdoubrud; zeros(4*fs,1)];
statesdoubrud = [delvpaddoubrud, alppaddoubrud, Qpaddoubrud, thetapaddoubrud];
inputsdoubrud = [depaddoubrud, drpaddoubrud];
ydoubrud = [statesdoubrud inputsdoubrud];
% pad the time signal
t4paddoubrud = 0:wsample:4-wsample;
t4paddoubrud = t4paddoubrud';
timepaddoubrud = [t4paddoubrud+timedoubrud(1)-t4paddoubrud(end)-wsample; timedoubrud; t4paddoubrud+timedoubrud(end)+wsample]; 
[yrsdoubrud,trsdoubrud] = resample(ydoubrud,timepaddoubrud,fs);

% plot to check resampling
figure, plot(trsdoubrud,yrsdoubrud,'x'), hold on, plot(timepaddoubrud,ydoubrud,'o')

% chop off padded data
%7621 is the last data point
yrsdoubrud(1:length(t4paddoubrud)+2,:)=[];
yrsdoubrud(end-length(t4paddoubrud)-1:end,:)=[];
trsdoubrud(1:length(t4paddoubrud)+2,:)=[];
trsdoubrud(end-length(t4paddoubrud)-1:end,:)=[];
trsdoubrud = trsdoubrud-trsdoubrud(1);

%% gather data for sys id

Q = ttest.Qrad_s(26075:30550);
P = ttest.Prad_s(26075:30550);
R = ttest.Rrad_s(26075:30550);
theta = ttest.pitchdeg(26075:30550)*pi/180-.44*pi/180;
phi = ttest.rolldeg(26075:30550)*pi/180;
alpha = ttest.alphadeg(26075:30550)*pi/180-.2*pi/180;
beta = ttest.betadeg(26075:30550)*pi/180;
vtrue = ttest.Vtruektas(26075:30550)*1.68781-172;
v=beta.*vtrue;
w=alpha.*vtrue;
U=sqrt(vtrue.^2-v.^2-w.^2);
de=ttest.elev1deg(26075:30550)*pi/180; % de 172 +-24�
de = de-1.029*pi/180;
da=(ttest.Railn1deg(26075:30550)/2+(-1*ttest.Lailn1deg(26075:30550)/2))*pi/180; %da 172 +-15
dr=ttest.rudr1deg(26075:30550)*pi/180; %dr 172 +-17
dt=ttest.thro1part(26075:30550);
% h=ttest.alt1ftmsl;
time=ttest.realtime(26075:30550);
timetot=ttest.totltime(26075:30550);
time = time - time(1);


% plot data
figure
subplot(5,1,1)
plot(time,beta*180/pi,'lineWidth',2)
xlabel('Time (s)')
ylabel({'Sideslip';'(deg)'})
xlim([0 71])
subplot(5,1,2)
plot(time,P*180/pi,'lineWidth',2)
xlabel('Time (s)')
ylabel({'Roll Rate';'(deg/s)'})
xlim([0 71])
subplot(5,1,3)
plot(time,R*180/pi,'lineWidth',2)
xlabel('Time (s)')
ylabel({'Yaw Rate';'(deg/s)'})
xlim([0 71])
subplot(5,1,4)
plot(time,phi*180/pi,'lineWidth',2)
xlabel('Time (s)')
ylabel({'Roll Angle';'(deg)'})
xlim([0 71])
subplot(5,1,5)
plot(time,da*180/pi,'lineWidth',2),hold on,plot(time,dr*180/pi,'lineWidth',2)
xlabel('Time (s)')
ylabel({'Surface';'Deflection';'(deg)'})
xlim([0 71])
legend('Aileron','Rudder')

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
delvpad = [ zeros(4*fs,1); beta; zeros(4*fs,1)];
alppad = [ zeros(4*fs,1); P; zeros(4*fs,1)];
Qpad = [ zeros(4*fs,1); R; zeros(4*fs,1)];
thetapad = [ zeros(4*fs,1); phi; zeros(4*fs,1)];
depad = [ zeros(4*fs,1); da; zeros(4*fs,1)];
drpad = [ zeros(4*fs,1); dr; zeros(4*fs,1)];
states = [delvpad, alppad, Qpad, thetapad];
inputs = [depad, drpad];
y = [states inputs];
% pad the time signal
t4pad = 0:wsample:4-wsample;
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
U = yrs(:,5:6);
Xdot = deriv(X,1/fs);

% 0 - longitudinal states
% 1 - lateral-directional states
% 2 - all states
modelType = 1;

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
% p0 = [-.2;-0.001;1;-.0001;3;-2;.01;.1;-2;.001;-.1;-.001;0;1;0;0;0;0;-5;-.5;.1;1;0;0];
% [y,p,crb,rr]=oe('oelatdir',p0,U,trs,x0,0,X);
% 
% Aoe = [p(1),p(2),p(3),p(4);
%    p(5),p(6),p(7),p(8)
%    p(9),p(10),p(11),p(12)
%    0,1,0,0];
% Boe = [p(17),p(18);p(19),p(20);p(21),p(22);p(23),p(24)];

Aoe = [-0.288072189242682,-0.0533162818169352,0.988558697807192,-0.176861075000975;3.41747826142298,-2.15044527818237,0.446623508340006,0.00182130920393708;-2.13677740613900,-0.00131889451261433,-0.242904067401692,-0.00511478037241558;0,1,0,0];
Boe = [0.0108677767040877,0.0384166872422263;-5.90214416386603,-0.604433132497319;0.0891949918802059,1.43142287754156;-0.00479994583830369,-0.0298750066489124];
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

%% compare ail doub
figure
pzmap(sysOE,sysID)

x0 = [yrsdoubail(1,1) yrsdoubail(1,2) yrsdoubail(1,3) yrsdoubail(1,4)];
lsimoutdoubOE = lsim(sysOE,yrsdoubail(:,5:6),trsdoubail,x0);
lsimoutdoubSINDY = lsim(sysID,yrsdoubail(:,5:6),trsdoubail,x0);
figure, 
subplot(5,1,1)
plot(trsdoubail,yrsdoubail(:,1)*180/pi,'lineWidth',2),hold on, plot(trsdoubail,lsimoutdoubOE(:,1)*180/pi,'lineWidth',2),plot(trsdoubail,lsimoutdoubSINDY(:,1)*180/pi,'lineWidth',2)
xlabel('Time (s)')
ylabel({'Sideslip';'(deg)'})
xlim([0 10])

legend('X-Plane data','Output Error Method','SINDY Method')
subplot(5,1,2)
plot(trsdoubail,yrsdoubail(:,2)*180/pi,'lineWidth',2),hold on, plot(trsdoubail,lsimoutdoubOE(:,2)*180/pi,'lineWidth',2),plot(trsdoubail,lsimoutdoubSINDY(:,2)*180/pi,'lineWidth',2)
xlabel('Time (s)')
ylabel({'Roll Rate';'(deg/s)'})
xlim([0 10])
subplot(5,1,3)
plot(trsdoubail,yrsdoubail(:,3)*180/pi,'lineWidth',2),hold on, plot(trsdoubail,lsimoutdoubOE(:,3)*180/pi,'lineWidth',2),plot(trsdoubail,lsimoutdoubSINDY(:,3)*180/pi,'lineWidth',2)
xlabel('Time (s)')
ylabel({'Yaw Rate';'(deg/s)'})
xlim([0 10])
subplot(5,1,4)
plot(trsdoubail,yrsdoubail(:,4)*180/pi,'lineWidth',2),hold on, plot(trsdoubail,lsimoutdoubOE(:,4)*180/pi,'lineWidth',2),plot(trsdoubail,lsimoutdoubSINDY(:,4)*180/pi,'lineWidth',2)
xlabel('Time (s)')
ylabel({'Roll Angle';'(deg)'})
xlim([0 10])
subplot(5,1,5)
plot(trsdoubail,yrsdoubail(:,5)*180/pi,'lineWidth',2),hold on, plot(trsdoubail,yrsdoubail(:,6)*180/pi,'lineWidth',2)
xlabel('Time (s)')
ylabel({'Surface';'Deflection';'(deg)'})
legend('Aileron','Rudder')
xlim([0 10])
%% compare Rud doub
% trsdoubrud
% yrsdoubrud

x0 = [yrsdoubrud(1,1) yrsdoubrud(1,2) yrsdoubrud(1,3) yrsdoubrud(1,4)];
lsimoutdoubOE = lsim(sysOE,yrsdoubrud(:,5:6),trsdoubrud,x0);
lsimoutdoubSINDY = lsim(sysID,yrsdoubrud(:,5:6),trsdoubrud,x0);
figure, 
subplot(5,1,1)
plot(trsdoubrud,yrsdoubrud(:,1)*180/pi,'lineWidth',2),hold on, plot(trsdoubrud,lsimoutdoubOE(:,1)*180/pi,'lineWidth',2),plot(trsdoubrud,lsimoutdoubSINDY(:,1)*180/pi,'lineWidth',2)
xlabel('Time (s)')
ylabel({'Sideslip';'(deg)'})
xlim([0 11.5])

legend('X-Plane data','Output Error Method','SINDY Method')
subplot(5,1,2)
plot(trsdoubrud,yrsdoubrud(:,2)*180/pi,'lineWidth',2),hold on, plot(trsdoubrud,lsimoutdoubOE(:,2)*180/pi,'lineWidth',2),plot(trsdoubrud,lsimoutdoubSINDY(:,2)*180/pi,'lineWidth',2)
xlabel('Time (s)')
ylabel({'Roll Rate';'(deg/s)'})
xlim([0 11.5])
subplot(5,1,3)
plot(trsdoubrud,yrsdoubrud(:,3)*180/pi,'lineWidth',2),hold on, plot(trsdoubrud,lsimoutdoubOE(:,3)*180/pi,'lineWidth',2),plot(trsdoubrud,lsimoutdoubSINDY(:,3)*180/pi,'lineWidth',2)
xlabel('Time (s)')
ylabel({'Yaw Rate';'(deg/s)'})
xlim([0 11.5])
subplot(5,1,4)
plot(trsdoubrud,yrsdoubrud(:,4)*180/pi,'lineWidth',2),hold on, plot(trsdoubrud,lsimoutdoubOE(:,4)*180/pi,'lineWidth',2),plot(trsdoubrud,lsimoutdoubSINDY(:,4)*180/pi,'lineWidth',2)
xlabel('Time (s)')
ylabel({'Roll Angle';'(deg)'})
xlim([0 11.5])
subplot(5,1,5)
plot(trsdoubrud,yrsdoubrud(:,5)*180/pi,'lineWidth',2),hold on, plot(trsdoubrud,yrsdoubrud(:,6)*180/pi,'lineWidth',2)
xlabel('Time (s)')
ylabel({'Surface';'Deflection';'(deg)'})
legend('Aileron','Rudder')
xlim([0 11.5])

%% Calculate Derivs

m = 13000/32.2;
u0 = 182;
theta0 = 0.05;
Ixx = 28000;
Izz = 47000;
Izx = 1300;
Ixz = Izx;
qbar = 34.3
rho = 0.0023;
S = 230;
b = 34;
A_sysID(1,2:4)=A_sysID(1,2:4)*-1;
A_sysID(2:4,1) = A_sysID(2:4,1)*-1;


[Coefs] = latdirStabDerivs(sysID.A,sysID.B,m, u0, theta0, Ixx, Izz, Izx, rho, S, b)
Coefs.Properties.VariableNames = {'Sindy'}





%% Calculate derivatives Roskam
qbar = 134.6;
S = 230;
m = 13000;
U1 = 182;
CD1 = 0.0279;
CL1 = 0.28;
cbar = 7;
Iyy = 18800;


clbeta = -0.173;
clp = -0.39;
clr = 0.45;
cybeta = -0.73;
cyp = 0;
cyr = 0.4;
cnbeta = 0.15;
cnp = -0.13;
cnr = -0.26;

clda = 0.149;
cldr = 0.014;
cyda = 0;
cydr = 0.14;
cnda = -0.05;
cndr = -0.074;

betadotbeta = qbar*S/(m*U1)*cybeta;
betadotp = 0.5*qbar*S*b/(m*U1^2)*cyp;
betadotr = 0.5*qbar*S*b/(m*U1^2)-1;
betadotphi = 32.17/U1;

ldden = Ixx*Izz-Ixz^2;

pdotbeta = (Izz*qbar*S*b*clbeta +Ixz*qbar*S*b*cnbeta)/ldden;
pdotp = (Izz*0.5*qbar*S*b^2/U1*clp+Ixz*0.5*qbar*S*b^2/U1*cnp)/ldden;
pdotr = (Izz*0.5*qbar*S*b^2/U1*clr+Ixz*0.5*qbar*S*b^2/U1*cnr)/ldden;
pdotphi = 0;

rdotbeta = (Ixx*qbar*S*b*cnbeta +Ixz*qbar*S*b*clbeta)/ldden;
rdotp = (Ixx*0.5*qbar*S*b^2/U1*cnp+Ixz*0.5*qbar*S*b^2/U1*clp)/ldden;
rdotr = (Ixx*0.5*qbar*S*b^2/U1*cnr+Ixz*0.5*qbar*S*b^2/U1*clr)/ldden;
rdotphi = 0;

AR = [betadotbeta betadotp betadotr betadotphi;
      pdotbeta pdotp pdotr pdotphi;
      rdotbeta rdotp rdotr rdotphi;
      0 1 0 0];

betadotda = qbar*S/(m*U1)*cyda;
betadotdr = qbar*S/(m*U1)*cydr;

pdotda = (Izz*qbar*S*b*clda +Ixz*qbar*S*b*cnda)/ldden;
pdotdr = (Izz*qbar*S*b*cldr +Ixz*qbar*S*b*cndr)/ldden;

rdotda = (Ixx*qbar*S*b*cnda +Ixz*qbar*S*b*clda)/ldden;
rdotdr = (Ixx*qbar*S*b*cndr +Ixz*qbar*S*b*cldr)/ldden;

BR = [betadotda betadotdr;
      pdotda pdotdr
      rdotda rdotdr
      0 0];

AR(1,2:4)=AR(1,2:4)*-1;
AR(2:4,1) = AR(2:4,1)*-1;

BR(2:4,1:2) = BR(2:4,1:2)*-1;

sysR = ss(AR,BR,C_sysID,D_sysID,'StateName',stateNames,'InputName',inputNames);
[CoefsR] = latdirStabDerivs(AR,BR,m, u0, theta0, Ixx, Izz, Izx, rho, S, b)
CoefsR.Properties.VariableNames = {'Roskam'}

%% Sindy
XuSindy = sysID.A(1,1);
XalpSindy = sysID.A(1,2);
XqSindy = sysID.A(1,3);
XdeSindy = sysID.B(1,1);
CDuSindy=((XuSindy*m*U1)/(-qbar*S))-2*CD1;
CDalpSindy=((XalpSindy*m)/(-qbar*S))+CL1;
CDqSindy = (XqSindy*m)/(-qbar*S);
CDdeSindy=((XdeSindy*m)/(-qbar*S));

ZuSindy = sysID.A(2,1);
ZalpSindy = sysID.A(2,2);
ZqSindy = sysID.A(2,3);
ZdeSindy = sysID.B(2,1);
CLuSindy=((ZuSindy*m*U1)/(-qbar*S))-2*CL1;
CLalpSindy=((ZalpSindy*m)/(-qbar*S))-CD1;
CLqSindy = 2*m*U1*ZqSindy/(-qbar*S*cbar);
CLdeSindy = ZdeSindy*m/(-qbar*S);

MuSindy = sysID.A(3,1);
MalpSindy = sysID.A(3,2);
MqSindy = sysID.A(3,3);
MdeSindy = sysID.B(3,1);

CmuSindy = (MuSindy*Iyy*U1)/(qbar*S*cbar)-(2*Cm1);
CmalpSindy = MalpSindy*Iyy/(qbar*S*cbar);
CmqSindy = (MqSindy*2*Iyy*U1)/(qbar*S*cbar^2);
CmdeSindy = MdeSindy*Iyy/(qbar*S*cbar);

%% OE
XuOE= sysOE.A(1,1);
XalpOE = sysOE.A(1,2);
XqOE = sysOE.A(1,3);
XdeOE = sysOE.B(1,1);
CDuOE=((XuOE*m*U1)/(-qbar*S))-2*CD1;
CDalpOE=((XalpOE*m)/(-qbar*S))+CL1;
CDqOE = (XqOE*m)/(-qbar*S);
CDdeOE=((XdeOE*m)/(-qbar*S));

ZuOE = sysOE.A(2,1);
ZalpOE = sysOE.A(2,2);
ZqOE = sysOE.A(2,3);
ZdeOE = sysOE.B(2,1);
CLuOE=((ZuOE*m*U1)/(-qbar*S))-2*CL1;
CLalpOE=((ZalpOE*m)/(-qbar*S))-CD1;
CLqOE = 2*m*U1*ZqOE/(-qbar*S*cbar);
CLdeOE = ZdeOE*m/(-qbar*S);

MuOE = sysOE.A(3,1);
MalpOE = sysOE.A(3,2);
MqOE = sysOE.A(3,3);
MdeOE = sysOE.B(3,1);

CmuOE = (MuOE*Iyy*U1)/(qbar*S*cbar)-(2*Cm1);
CmalpOE = MalpOE*Iyy/(qbar*S*cbar);
CmqOE = (MqOE*2*Iyy*U1)/(qbar*S*cbar^2);
CmdeOE = MdeOE*Iyy/(qbar*S*cbar);

[CoefsOE] = latdirStabDerivs(sysOE.A,sysOE.B,m, u0, theta0, Ixx, Izz, Izx, rho, S, b)
CoefsOE.Properties.VariableNames = {'Output_Error'}

LonStab = [CoefsR Coefs CoefsOE]
writetable(LonStab,'LearApproachLatDirStabDerivs.xlsx','WriteRowNames',true)

% Xu=A_sysID(1,1);
% CDu=((Xu*m*U1)/(-qbar*S))-2*CD1;
% 
% Xalp=A_sysID(1,2)*U1;
% CDalp=((Xalp*m)/(-qbar*S))+CL1;
% 
% Xde=B_sysID(1,1);
% CDde=((Xde*m)/(-qbar*S));
% 
% Zu=A_sysID(2,1);
% CLu=((Zu*m*U1)/(-qbar*S))-2*CL1;
% 
% Zalp=A_sysID(2,2)*U1;
% CLalp=((Zalp*m)/(-qbar*S))-CD1;
% 
% Zq=A_sysID(2,3);
% CLq=2*m*U1*Zq/(-qbar*S*cbar);
% 
% Zde=B_sysID(2,1);
% CLde=Zde*m/(-qbar*S);

%% Plots 


figure, 
pzmap(sysR,sysID,sysOE)
legend('Roskam','Sindy','Output-Error')
grid on

%% bodes

% sideslip aileron
figure
fig = gcf
fig.Units = 'inches' 
fig.Position = [3 2 8 5]
fig.Color = [1 1 1]
lonbodedu = bodeplot(sysID,sysOE,sysR)

lonbodedu.OutputName = {'Sideslip','Roll Rate','Yaw Rate','Roll Angle'};
lonbodedu.OutputVisible = {'on','off','off','off'}
lonbodedu.InputVisible = {'on','off'}
xlim([0.1 10])
title('')
grid on
legend('SINDY','Output-Error','Roskam','Position',[0.75 0.75 0.1 0.1])
lonbodedu.AxesGrid.RowLabel = ''


p = getoptions(lonbodedu); 
p.PhaseWrapping = 'on' 
setoptions(lonbodedu,p)
lonbodedu.AxesGrid.YLabel = {'Phase','Magnitude'}
lonbodedu.AxesGrid.YUnits = {'Deg','DB'}

% p ail
figure
fig = gcf
fig.Units = 'inches' 
fig.Position = [3 2 8 5]
fig.Color = [1 1 1]
lonbodedu = bodeplot(sysID,sysOE,sysR)

lonbodedu.OutputName = {'\Deltau','Angle of Attack','Pitch Rate','Pitch Angle'};
lonbodedu.OutputVisible = {'off','on','off','off'}
lonbodedu.InputVisible = {'on','off'}
xlim([0.1 10])

title('')
grid on
legend('SINDY','Output-Error','Roskam','Position',[0.6 0.52 0.1 0.1])
lonbodedu.AxesGrid.RowLabel = ''

lonbodedu.AxesGrid.YLabel = {'Phase','Magnitude'}
lonbodedu.AxesGrid.YUnits = {'Deg','DB'}

% r ail
figure
fig = gcf
fig.Units = 'inches' 
fig.Position = [3 2 8 5]
fig.Color = [1 1 1]
lonbodedu = bodeplot(sysID,sysOE,sysR)

lonbodedu.OutputName = {'\Deltau','Angle of Attack','Pitch Rate','Pitch Angle'};
lonbodedu.OutputVisible = {'off','off','on','off'}
lonbodedu.InputVisible = {'on','off'}
xlim([0.1 10])
title('')
grid on
legend('SINDY','Output-Error','Roskam','Position',[0.77 0.52 0.1 0.1])
lonbodedu.AxesGrid.RowLabel = ''


p = getoptions(lonbodedu); 
p.PhaseWrapping = 'on' 
setoptions(lonbodedu,p)
lonbodedu.AxesGrid.YLabel = {'Phase','Magnitude'}
lonbodedu.AxesGrid.YUnits = {'Deg','DB'}

%% rudder

% sideslip aileron
figure
fig = gcf
fig.Units = 'inches' 
fig.Position = [3 2 8 5]
fig.Color = [1 1 1]
lonbodedu = bodeplot(sysID,sysOE,sysR)

lonbodedu.OutputName = {'Sideslip','Roll Rate','Yaw Rate','Roll Angle'};
lonbodedu.OutputVisible = {'on','off','off','off'}
lonbodedu.InputVisible = {'off','on'}
xlim([0.1 10])
title('')
grid on
legend('SINDY','Output-Error','Roskam','Position',[0.75 0.75 0.1 0.1])
lonbodedu.AxesGrid.RowLabel = ''


p = getoptions(lonbodedu); 
p.PhaseWrapping = 'on' 
setoptions(lonbodedu,p)
lonbodedu.AxesGrid.YLabel = {'Phase','Magnitude'}
lonbodedu.AxesGrid.YUnits = {'Deg','DB'}

% p ail
figure
fig = gcf
fig.Units = 'inches' 
fig.Position = [3 2 8 5]
fig.Color = [1 1 1]
lonbodedu = bodeplot(sysID,sysOE,sysR)

lonbodedu.OutputName = {'\Deltau','Angle of Attack','Pitch Rate','Pitch Angle'};
lonbodedu.OutputVisible = {'off','on','off','off'}
lonbodedu.InputVisible = {'off','on'}
xlim([0.1 10])
title('')
grid on
legend('SINDY','Output-Error','Roskam','Position',[0.6 0.52 0.1 0.1])
lonbodedu.AxesGrid.RowLabel = ''

lonbodedu.AxesGrid.YLabel = {'Phase','Magnitude'}
lonbodedu.AxesGrid.YUnits = {'Deg','DB'}

% r ail
figure
fig = gcf
fig.Units = 'inches' 
fig.Position = [3 2 8 5]
fig.Color = [1 1 1]
lonbodedu = bodeplot(sysID,sysOE,sysR)

lonbodedu.OutputName = {'\Deltau','Angle of Attack','Pitch Rate','Pitch Angle'};
lonbodedu.OutputVisible = {'off','off','on','off'}
lonbodedu.InputVisible = {'off','on'}
xlim([0.1 10])
title('')
grid on
legend('SINDY','Output-Error','Roskam','Position',[0.77 0.52 0.1 0.1])
lonbodedu.AxesGrid.RowLabel = ''


p = getoptions(lonbodedu); 
p.PhaseWrapping = 'on' 
setoptions(lonbodedu,p)
lonbodedu.AxesGrid.YLabel = {'Phase','Magnitude'}
lonbodedu.AxesGrid.YUnits = {'Deg','DB'}




















