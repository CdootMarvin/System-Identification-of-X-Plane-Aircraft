addpath('C:\Users\Chris\Documents\THESIShome\Sparse Identification\sparsedynamics\utils')
addpath('C:\Users\Chris\Documents\THESIShome\SIDPAC_ver_4.1\SIDPAC')

clear all
close all
clc

%% Linear Models
% Lon


  Alon = [-0.0008223      -13.96     -0.9411      -31.96;
  -0.0001764      -1.653      0.9912           0
     0.0006181      -25.05      -1.129     -0.0181
           0           0      0.9981           0];
 
  Blon = [2.939
-0.1261
-29.58
0];
 
  Clon = eye(4);
  Dlon = zeros(4,1);
  
Klon = [0    0.4198   -0.1663   -0.0210];
  
  Acl = Alon-Blon*Klon;
  
  sysCLlon = ss(Acl,Blon,Clon,Dlon);
  syslon = ss(Alon,Blon,Clon,Dlon);
%====================================================
% Preprocess the text file so we can read it in with readtable().
fInput = fopen('LearCruiseOpenVClose.txt', 'rt');
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

Q = ttest.Qrad_s;
P = ttest.Prad_s;
R = ttest.Rrad_s;
theta = ttest.pitchdeg*pi/180-.44*pi/180;
phi = ttest.rolldeg*pi/180;
alpha = ttest.alphadeg*pi/180-.2*pi/180;
beta = ttest.betadeg*pi/180;
vtrue = ttest.Vtruektas*1.68781;
v=beta.*vtrue;
w=alpha.*vtrue;
U=sqrt(vtrue.^2-v.^2-w.^2);
de=ttest.elev1deg*pi/180; % de 172 +-24°
%de = de-1.029*pi/180;
da=(ttest.Railn1deg/2+(-1*ttest.Lailn1deg/2))*pi/180; %da 172 +-15
dr=ttest.rudr1deg*pi/180; %dr 172 +-17
dt=ttest.thro1part;
hdg = ttest.hdingtrue;
% h=ttest.alt1ftmsl;
time=ttest.realtime;
timetot=ttest.totltime;
time = time - time(1);


% plot data lat dir
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

% plot data lon
figure
subplot(5,1,1)
plot(time,vtrue)
xlabel('time (s)')
ylabel('vtas (ft/s)')
subplot(5,1,2)
plot(time,alpha*180/pi)
xlabel('time (s)')
ylabel('Angle of Attack (deg)')
subplot(5,1,3)
plot(time,Q*180/pi)
xlabel('time (s)')
ylabel('Pitch Rate (deg/s)')
subplot(5,1,4)
plot(time,theta*180/pi)
xlabel('time (s)')
ylabel('Pitch Angle (deg)')
subplot(5,1,5)
plot(time,de*180/pi)
xlabel('time (s)')
ylabel('Surface Deflection (deg)')
legend('elevator')


%% Roll
close all
% Bare A/C
% time = 620 to 634 = 22838 tot 23280
% plot data lat dir
% Closed Loop
% time = 476 to 488 = 17732 to 18109
timerollbac = time(22838:23280);
timerollbac = timerollbac-timerollbac(1);
betarollbac = beta(22838:23280);
prollbac = P(22838:23280);
rrollbac = R(22838:23280);
phirollbac = phi(22838:23280);
darollbac = da(22838:23280);
drrollbac = dr(22838:23280);

timerollcl = time(17732:18109);
timerollcl = timerollcl-timerollcl(1);
betarollcl = beta(17732:18109);
prollcl = P(17732:18109);
rrollcl = R(17732:18109);
phirollcl = phi(17732:18109);
darollcl = da(17732:18109);
drrollcl = dr(17732:18109);

figure
subplot(5,1,1)
plot(timerollbac,betarollbac*180/pi,timerollcl,betarollcl*180/pi)
xlabel('time (s)')
ylabel('beta (deg)')
subplot(5,1,2)
plot(timerollbac,prollbac*180/pi,timerollcl,prollcl*180/pi)
xlabel('time (s)')
ylabel('Roll Rate (deg/s)')
subplot(5,1,3)
plot(timerollbac,rrollbac*180/pi,timerollcl,rrollcl*180/pi)
xlabel('time (s)')
ylabel('Yaw Rate (deg/s)')
subplot(5,1,4)
plot(timerollbac,phirollbac*180/pi,timerollcl,phirollcl*180/pi)
xlabel('time (s)')
ylabel('Roll Angle (deg)')
subplot(5,1,5)
plot(timerollbac,darollbac*180/pi)%,'Color',"#0072BD", 'LineStyle','--')
hold on,
plot(timerollbac,drrollbac*180/pi)%,'Color',"#0072BD", 'LineStyle','-')
plot(timerollcl,darollcl*180/pi)%,'Color',"#D95319", 'LineStyle','--')
plot(timerollcl,drrollcl*180/pi)%,'Color',"#D95319", 'LineStyle','-')
xlabel('time (s)')
ylabel('Surface Deflection (deg)')
legend('Aileron','Rudder')



%% Yaw
% Bare A/C
% time = 648 to 660 = 23821 to 24202
% Closed Loop
% time = 508 to 520 = 18875 to 19252

timeyawbac = time(23827:24202);
timeyawbac = timeyawbac-timeyawbac(1);
betayawbac = beta(23827:24202);
pyawbac = P(23827:24202);
ryawbac = R(23827:24202);
phiyawbac = phi(23827:24202);
dayawbac = da(23827:24202);
dryawbac = dr(23827:24202);

timeyawcl = time(18875:19252);
timeyawcl = timeyawcl-timeyawcl(1);
betayawcl = beta(18875:19252);
pyawcl = P(18875:19252);
ryawcl = R(18875:19252);
phiyawcl = phi(18875:19252);
dayawcl = da(18875:19252);
dryawcl = dr(18875:19252);

figure
subplot(5,1,1)
plot(timeyawbac,betayawbac*180/pi,timeyawcl,betayawcl*180/pi)
xlabel('time (s)')
ylabel('beta (deg)')
subplot(5,1,2)
plot(timeyawbac,pyawbac*180/pi,timeyawcl,pyawcl*180/pi)
xlabel('time (s)')
ylabel('Roll Rate (deg/s)')
subplot(5,1,3)
plot(timeyawbac,ryawbac*180/pi,timeyawcl,ryawcl*180/pi)
xlabel('time (s)')
ylabel('Yaw Rate (deg/s)')
subplot(5,1,4)
plot(timeyawbac,phiyawbac*180/pi,timeyawcl,phiyawcl*180/pi)
xlabel('time (s)')
ylabel('Roll Angle (deg)')
subplot(5,1,5)
plot(timeyawbac,dayawbac*180/pi)%,'Color',"#0072BD", 'LineStyle','--')
hold on,
plot(timeyawbac,dryawbac*180/pi)%,'Color',"#0072BD", 'LineStyle','-')
plot(timeyawcl,dayawcl*180/pi)%,'Color',"#D95319", 'LineStyle','--')
plot(timeyawcl,dryawcl*180/pi)%,'Color',"#D95319", 'LineStyle','-')
xlabel('time (s)')
ylabel('Surface Deflection (deg)')
legend('Aileron','Rudder')

%% Pitch
% Bare A/C
% time = 566 to 585 = 7336 to 7651
% Close Loop
% time 410 to 425 = 9855 to 10169
timelonbac = time(20547:20988);
timelonbac = timelonbac-timelonbac(1);
vlonbac = vtrue(20547:20988)/1.68781;
aoalonbac = alpha(20547:20988);
qlonbac = Q(20547:20988);
thetalonbac = theta(20547:20988);
delonbac = de(20547:20988);

timeloncl = time(14757:15266);
timeloncl = timeloncl-timeloncl(1);
vloncl = vtrue(14757:15266)/1.68781;
aoaloncl = alpha(14757:15266);
qloncl = Q(14757:15266);
thetaloncl = theta(14757:15266);
deloncl = de(14757:15266);

for ind = 1:length(timelonbac)-1
    delt(ind) = timelonbac(ind+1)-timelonbac(ind);
    ind = ind+1;
end
wsample = mean(delt);
fs = round(1/wsample);

depad = [ zeros(4*fs,1); delonbac; zeros(4*fs,1)];
t4paddoub = 0:wsample:4-2*wsample;
t4paddoub = t4paddoub';
timepaddoub = [t4paddoub+timelonbac(1)-t4paddoub(end)-wsample; timelonbac; t4paddoub+timelonbac(end)+wsample]; 
[ders,trs] = resample(depad,timepaddoub,fs);

ders(1:length(t4paddoub),:)=[];
ders(end-length(t4paddoub):end,:)=[];
trs(1:length(t4paddoub),:)=[];
trs(end-length(t4paddoub):end,:)=[];


figure, plot(trs,ders,'x',timeloncl,deloncl,'o')
ylon = lsim(sysCLlon, ders, trs, [0, aoaloncl(1)-0.9*pi/180, qloncl(1), thetaloncl(1)]);


figure, plot(trs,ylon(:,2)*180/pi+0.9,timelonbac,aoalonbac*180/pi,timeloncl,aoaloncl*180/pi)
figure, plot(trs,ylon(:,3)*180/pi,timelonbac,qlonbac*180/pi,timeloncl,qloncl*180/pi)
figure, plot(trs, ylon(:,1)/1.68781+vloncl(1),timelonbac,vlonbac,timeloncl,vloncl)
figure, plot(trs, ylon(:,4)*180/pi,timelonbac,thetalonbac*180/pi,timeloncl,thetaloncl*180/pi)

figure
subplot(5,1,1)
plot(timelonbac,vlonbac,timeloncl,vloncl,trs, ylon(:,1)/1.68781+vloncl(1))
xlabel('time (s)')
ylabel('vtas (Kts)')
subplot(5,1,2)
plot(timelonbac,aoalonbac*180/pi,timeloncl,aoaloncl*180/pi,trs,ylon(:,2)*180/pi+0.9)
xlabel('time (s)')
ylabel('Angle of Attack (deg)')
subplot(5,1,3)
plot(timelonbac,qlonbac*180/pi,timeloncl,qloncl*180/pi,trs,ylon(:,3)*180/pi)
xlabel('time (s)')
ylabel('Pitch Rate (deg/s)')
subplot(5,1,4)
plot(timelonbac,thetalonbac*180/pi,timeloncl,thetaloncl*180/pi,trs,ylon(:,4)*180/pi)
xlabel('time (s)')
ylabel('Pitch Angle (deg)')
subplot(5,1,5)
plot(timelonbac,delonbac*180/pi,timeloncl,deloncl*180/pi)
xlabel('time (s)')
ylabel('Surface Deflection (deg)')
legend('elevator')


