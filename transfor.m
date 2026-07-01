clear
clc
Ta=22;%ambient temperchare;
T_tor=55;%top oil temprature rise at rated load
P_rated=500;%rated load(MVA)
P_actual =140;%Actual load(MVA)
T_G=50;%winding hottest-spot temperature rise over top-oil temperature.
R_L=4.1;%Loss ratio
K=(P_actual/P_rated);
T_HS=Ta+T_G+T_tor*(((((K^2)*R_L)+1)/(R_L+1)))^.92;
%%%%%%%%%%%%%%%%%%%%%%%%
F=.0008;%input('please Enetr furfural Value ? ');
if F>=.5
    F=4;
end
% The activation energy
Ea = 111*10^3; %J/mol
% Process Constant
A = 2*10^8; %h-1
% Gas Costant
Rg = 8.314;
%Temperature (As a function of time. Assumed constant)
T = T_HS+273; %K, 98 deg C
% Reaction rate
k = A*exp(-Ea/(Rg*T)); %h-1
k = 24*k; %d-1
% Original DP value (Degree of Polymerization)
DP0 = (1.51-log10(F))/0.00355;
% Threshold DP value below which transformer ceases to function
 DPc =250;
%% Create the ideal transformer
% Start of history
t0 = 1;
% End of history
tm = 45*365; % Arbitrarily chosen as 15 years since average transformer life with above parameters was aobserved to about 10-12 years
% Generate the curve
DP(t0) = DP0;
cur=0;
for t=t0:tm
DP(t) = DP(t0)/(1+DP(t0)*k*(t-t0));


if DP(t)<=DPc
    cur=cur+1;
    endlife(cur)=t/365;
end
end
endlife1=(endlife(1))
%%
cur1=0;
%%%%%%%%%%%%%%%%% Ideal transform %%%%%%%%%%%%%%%
sp0 = -2; % shape parameter
sig0 = 0.02;
mu0 =0;
% Generate deviations for each time-step
for i=t0:tm
Dev0(i) = gevrnd(sp0, sig0, mu0);
end
% figure,
% hist(Dev0);
% Generate the curve by introducing deviations to the original curve
DPa(t0) = DP0;
for t=t0+1:tm
DPa(t) = DPa(t-1)/(1+DPa(t-1)*k*1);
DPa(t) = max( DPa(t) + min(Dev0(t),0), 0);
if DPa(t) < DPc && DPa(t-1) >= DPc
disp(sprintf('Non-ideal transformer would last for %d days or %d years',t, t/365));
TL=t;
end
if DPa(t)<=DPc
    cur1=cur1+1;
    endlife(cur1)=t/365;
end
end
endlife2=(endlife(1))

% plot((t0:tm)./365, DPa, (t0:tm)./365, DP);
% plot((t0:tm)./365, DPa);
% title('Ideal and Actual Transformer History Curve');
% xlabel('Number of days');
% ylabel('DP value');
%%%%%%%%%%%%%%%%  Compute the Failure %%%%%%%%%%%%%%%%%%%%%%
%% Compute the probability that transformer will fail in next m days
% Prediction days
m = 10;
% Today
%today = TL-5m:TL;
Tw = 15;
for today=TL-Tw:TL
%% Compute the distribution of deviations based on history
for t=t0+1:today
ind = find(DP<DPa(t-1),1);
Dev(t-1) = DPa(t) - DP(ind);
end
%figure, hist(Dev);
paramhat = gevfit(Dev);
sp = paramhat(1); sig = paramhat(2); mu = paramhat(3);
% disp(sprintf('Curve fit successful ! sp = %d scale = %d mu = %d ',sp, sig, mu));
%% Rgun the simulation for next m days
% Number of iterations
numit = 10000;
itcount = 0;
for it=1:numit
currentDP = DPa(today);
for tn=1:m
ind = find(DP<currentDP,1);
currentDP = max( DP(ind) + min(gevrnd(sp,sig,mu),0),0);
if currentDP < DPc
itcount = itcount + 1;
break;
end
end
end
pf = itcount/numit;
%disp(sprintf('Probability of failure on day %d is %d', today, pf));
pfv(today) = pf;
clc
end
clc
figure, 
subplot(2,1,1)
plot((t0:tm)./365, DPa);
title('Ideal and Actual Transformer History Curve');
xlabel('Number of Years');
ylabel('DP value');
subplot(2,1,2)
plot((TL-Tw:TL)./365,pfv(TL-Tw:TL));
title('Failure Probability Across Transformer Lifetime');
xlabel('Number of Years');
ylabel('Probability of failure');

figure
Ta=45;%ambient temperchare;
T_tor=45;%top oil temprature rise at rated load
P_rated=500;%rated load(MVA)
P_actual =50:10:500*1.3;%Actual load(MVA)
T_G=50;%winding hottest-spot temperature rise over top-oil temperature.
R_L=4.1;%Loss ratio
K=(P_actual/P_rated);
T_HS=Ta+T_G+T_tor.*(((((K.^2)*R_L)+1)./(R_L+1))).^0.92;
plot(P_actual,T_HS)
xlabel('Actual load(MVA)');
ylabel('Temperature (Deg)');
