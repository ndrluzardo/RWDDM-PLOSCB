% Extinction with CS time change
% Phase 1: simple acquisition
% Phase 2: extinction with a longer CS
% phase 3: extinction with a shorter CS
% phase 4: extinction with the same duration CS

clear
close all

%% Phase 1 - Acquisition

cycle_num=150;

CS_dur=20;

%---parameters
h=0.01;
tau_x=50;
alpha_t=0.09;
alpha_E=0.06;
mu=1;
sigma=0.35;
m=0.16;
H=30;
%---

%---structure arrays
CS.name='CS';
%-A
CS.A.ph1=zeros(cycle_num,1);
CS.A.ph2=[];
CS.A.ph1(1)=1*10^(-3);
%-V
CS.V.ph1=zeros(cycle_num,1);
CS.V.ph2=[];
%---

%---DDM constants
N=normrnd(0,1,ceil(CS_dur/h*cycle_num),1); % noise for CS A
%---

counterDDM=0;
cycle_length=round(CS_dur/h);

for trial=1:cycle_num
    
    %--initialize values for timer and CS
    P=zeros(1,cycle_length);
    x=zeros(1,cycle_length);
    %--
    
    for t=1:cycle_length
        
        counterDDM=counterDDM+1; % update counter for random process in DDM
        
        A=CS.A.ph1(trial);
        V=CS.V.ph1(trial);
        
        % min will take the minimum value: either DDM result or 3. This
        % caps the value of integrator at 3.
        P(t+1)=min(DDM( P(t), A, h, m, N(counterDDM) ), 3);
        
        % max ensures the minimum value the accumulator can reach is
        % 1*10^-6. This avoids division by zero later.
        P(t+1)=max(P(t+1), 1*10^(-3));
        
        %---Element (RBFs)
        x(t)=CStrace(P(t+1),mu,sigma,tau_x,1,x(t),h);
        %---
        
    end
    
    %---Slope Correction
    A=A+A*alpha_t*(1-P(t+1))/P(t+1); % realistic correction rule, never fully converges. Only updates in rewarded trials.
    %---
    
    %---V update
    V=RW(V,alpha_E,x(t),H,A,P(t+1));
    %---
    
    CS.A.ph1(trial+1)=A;
    CS.V.ph1(trial+1)=V;
    
end

plot(CS.V.ph1,'LineWidth',5)
PlotProperties


%% Phase 2 - Extinction with longer CS

%---Number of cycles (trials)
cycle_num=120;
%---

%---CS duration
CS_dur=40;
%---

%---structure arrays
%-A
CS.A.ph2=zeros(cycle_num,1);
CS.A.ph2(1)=CS.A.ph1(end);
%-V
CS.V.ph2=zeros(cycle_num,1);
CS.V.ph2(1)=CS.V.ph1(end);
%-CR
CS.CR.ph2=zeros(cycle_num,cycle_length);
%---

%---DDM constants
N=normrnd(0,1,ceil(CS_dur/h*cycle_num),1); % noise for CS A
%---

counterDDM=0;
cycle_length=round(CS_dur/h);

for trial=1:cycle_num
    
    %--initialize values for timer and CS
    P=zeros(1,cycle_length);
    x=zeros(1,cycle_length);
    %--
    
    for t=1:cycle_length
        
        counterDDM=counterDDM+1; % update counter for random process in DDM
        
        A=CS.A.ph2(trial);
        V=CS.V.ph2(trial);
        
        % min will take the minimum value: either DDM result or 3. This
        % caps the value of integrator at 3.
        P(t+1)=min(DDM( P(t), A, h, m, N(counterDDM) ), 3);
        
        % max ensures the minimum value the accumulator can reach is
        % 1*10^-6. This avoids division by zero later.
        P(t+1)=max(P(t+1), 1*10^(-3));
        
        %---Element (RBFs)
        x(t)=CStrace(P(t+1),mu,sigma,tau_x,1,x(t),h);
        %---
        
        %---CR
        CS.CR.ph2(trial,t)=x(t)*V;
        %---
        
    end
    
    %---Slope Correction
    A=A+A*alpha_t*(1-P(t+1))/P(t+1); % realistic correction rule, never fully converges. Only updates in rewarded trials.
    %---
    
    %---V update
    V=RW(V,alpha_E,x(t),0,A,P(t+1));
    %---
    
    CS.A.ph2(trial+1)=A;
    CS.V.ph2(trial+1)=V;
    
end

early_mean=mean(CS.CR.ph2(1:5,:));
middle_mean=mean(CS.CR.ph2(15:20,:));
late_mean=mean(CS.CR.ph2(40:45,:));
x_vector=h:h:length(early_mean)*h;

figure
plot(x_vector,early_mean,x_vector,middle_mean,x_vector,late_mean,'LineWidth',6)
%line([20 20], [0 max(early_mean)+0.5],'color','k','LineWidth',2)
legend({'early','middle','late'},'Location','northwest','Box','off')
title('model')
xlabel('time (sec)')
ylabel('response strength')
PlotProperties

figure
plot((1./CS.A.ph2)/1000,'LineWidth',6)
title('model')
xlabel('extinction trial')
ylabel('time estimate 1/A')
PlotProperties
%% Phase 3 - Extinction with shorter CS

%---Number of cycles (trials)
cycle_num=120;
%---

%---CS duration
CS_dur=10;
%---

%---Cycle length
cycle_length=round(CS_dur/h);
%---

%---structure arrays
%-A
CS.A.ph3=zeros(cycle_num,1);
CS.A.ph3(1)=CS.A.ph1(end);
%-V
CS.V.ph3=zeros(cycle_num,1);
CS.V.ph3(1)=CS.V.ph1(end);
%-CR
CS.CR.ph3=zeros(cycle_num,cycle_length);
%---

%---DDM constants
N=normrnd(0,1,ceil(CS_dur/h*cycle_num),1); % noise for CS A
%---

counterDDM=0;


for trial=1:cycle_num
    
    %--initialize values for timer and CS
    P=zeros(1,cycle_length);
    x=zeros(1,cycle_length);
    %--
    
    for t=1:cycle_length
        
        counterDDM=counterDDM+1; % update counter for random process in DDM
        
        A=CS.A.ph3(trial);
        V=CS.V.ph3(trial);
        
        % min will take the minimum value: either DDM result or 3. This
        % caps the value of integrator at 3.
        P(t+1)=min(DDM( P(t), A, h, m, N(counterDDM) ), 3);
        
        % max ensures the minimum value the accumulator can reach is
        % 1*10^-6. This avoids division by zero later.
        P(t+1)=max(P(t+1), 1*10^(-3));
        
        %---Element (RBFs)
        x(t)=CStrace(P(t+1),mu,sigma,tau_x,1,x(t),h);
        %---
        
        %---CR
        CS.CR.ph3(trial,t)=x(t)*V;
        %---
        
    end
    
    %---Slope Correction
    A=A+A*alpha_t*(1-P(t+1))/P(t+1); % realistic correction rule, never fully converges. Only updates in rewarded trials.
    %---
    
    %---V update
    V=RW(V,alpha_E,x(t),0,A,P(t+1));
    %---
    
    CS.A.ph3(trial+1)=A;
    CS.V.ph3(trial+1)=V;
    
end
early_mean2=mean(CS.CR.ph3(1:5,:));
middle_mean2=mean(CS.CR.ph3(15:20,:));
late_mean2=mean(CS.CR.ph3(40:45,:));
x_vector2=h:h:length(early_mean2)*h;

figure
plot(x_vector2,early_mean2,x_vector2,middle_mean2,x_vector2,late_mean2,'LineWidth',6)
%line([20 20], [0 max(early_mean)+0.5],'color','k','LineWidth',2)
legend({'early','middle','late'},'Location','northwest','Box','off')
xlabel('time (sec)')
ylabel('response strength')
PlotProperties

figure
plot((1./CS.A.ph3)/1000,'LineWidth',6)
xlabel('extinction trial')
ylabel('time estimate 1/A')
PlotProperties

%% Phase 4 - Extinction with same duration CS

%---Number of cycles (trials)
cycle_num=120;
%---

%---CS duration
CS_dur=20;
%---

%---Cycle length
cycle_length=round(CS_dur/h);
%---

%---structure arrays
%-A
CS.A.ph4=zeros(cycle_num,1);
CS.A.ph4(1)=CS.A.ph1(end);
%-V
CS.V.ph4=zeros(cycle_num,1);
CS.V.ph4(1)=CS.V.ph1(end);
%-CR
CS.CR.ph4=zeros(cycle_num,cycle_length);
%---

%---DDM constants
N=normrnd(0,1,ceil(CS_dur/h*cycle_num),1); % noise for CS A
%---

counterDDM=0;


for trial=1:cycle_num
    
    %--initialize values for timer and CS
    P=zeros(1,cycle_length);
    x=zeros(1,cycle_length);
    %--
    
    for t=1:cycle_length
        
        counterDDM=counterDDM+1; % update counter for random process in DDM
        
        A=CS.A.ph4(trial);
        V=CS.V.ph4(trial);
        
        % min will take the minimum value: either DDM result or 3. This
        % caps the value of integrator at 3.
        P(t+1)=min(DDM( P(t), A, h, m, N(counterDDM) ), 3);
        
        % max ensures the minimum value the accumulator can reach is
        % 1*10^-6. This avoids division by zero later.
        P(t+1)=max(P(t+1), 1*10^(-3));
        
        %---Element (RBFs)
        x(t)=CStrace(P(t+1),mu,sigma,tau_x,1,x(t),h);
        %---
        
        %---CR
        CS.CR.ph4(trial,t)=x(t)*V;
        %---
        
    end
    
    %---Slope Correction
    A=A+A*alpha_t*(1-P(t+1))/P(t+1); % realistic correction rule, never fully converges. Only updates in rewarded trials.
    %---
    
    %---V update
    V=RW(V,alpha_E,x(t),0,A,P(t+1));
    %---
    
    CS.A.ph4(trial+1)=A;
    CS.V.ph4(trial+1)=V;
    
end

plot(1:length(CS.V.ph2),CS.V.ph2,1:length(CS.V.ph3),CS.V.ph3,1:length(CS.V.ph4),CS.V.ph4,'LineWidth',5)
legend({'CS 40','CS 10','CS 20'},'Box','off')
title('model')
xlabel('extinction trial')
ylabel('associative strength')
PlotProperties

%% first 10-sec CR extinction analysis

% take only the first 10 seconds from each CS
CRlong=sum(CS.CR.ph2(:,1:round(10/h)),2);
CRshort=sum(CS.CR.ph3(:,1:round(10/h)),2);
CRsame=sum(CS.CR.ph4(:,1:round(10/h)),2);

%---Sum the CR over all time points, then apply a two-trial mean filter
A=zeros(120,60);
A(1:2,1)=1;
for k=2:60
    A(:,k)=circshift(A(:,k-1),2);
end
CRlong=(CRlong'*A)./200;
CRshort=(CRshort'*A)./200;
CRsame=(CRsame'*A)./200;
%---
figure
plot(1:length(CRlong),CRlong,1:length(CRshort),CRshort,1:length(CRsame),CRsame,'LineWidth',6)
legend({'CS 40','CS 10','CS 20'},'Box','off')
title('model')
xlabel('2 trial blocks')
ylabel('resp. strenght (arbitrary units)')
PlotProperties

save('ExtinctionTimeChangeSim','x_vector','x_vector2','early_mean','early_mean2','middle_mean','middle_mean2','late_mean','late_mean2','CS','CRlong','CRshort','CRsame')