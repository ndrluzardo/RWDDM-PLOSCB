% temporal averaging

clear
close all

%---Parameters
h=0.01;
tau_x=50;
alpha_t=0.2;
alpha_E=0.1;
mu=1;
sigma=0.35;
m=0.2;
H=30;
%---

CS_L_dur=20; % (ms)
CS_S_dur=10; % (ms)
PCS_S_dur=3*CS_S_dur; % peak CS short duration
PCS_L_dur=3*CS_L_dur; % peak CS long duration
PCS_C_dur=3.5*CS_L_dur; % peak CS compound duration

%--Create a random vector with trial durations
CS_dur=[CS_L_dur CS_S_dur]; % (ms)
CS_dur=repmat(CS_dur,1,700);
CS_dur=[CS_dur PCS_S_dur*ones(1,round((0.25/0.75)*length(CS_dur)/3)) PCS_L_dur*ones(1,round((0.25/0.75)*length(CS_dur)/3)) PCS_C_dur*ones(1,round((0.25/0.75)*length(CS_dur)/3))];
CS_dur=CS_dur(randperm(length(CS_dur)));
%---

cycle_num=length(CS_dur);

%---Create structure arrays to carry CS info
CS(1).name='short';
CS(2).name='long';
CS(1).rate=NaN(cycle_num,1);
CS(1).rate(1)=1*10^(-3);
CS(2).rate=NaN(cycle_num,1);
CS(2).rate(1)=1*10^(-3);
CS(1).V=zeros(cycle_num,1);
CS(2).V=zeros(cycle_num,1);
%---


%---DDM constants
N=normrnd(0,1,ceil(PCS_C_dur/h)*cycle_num,1); % noise for CS A
%---

%---Initialize counters
counterDDM=0;
counterCS_L=1;
counterCS_S=1;
%---

%---CR
CR=NaN(cycle_num,round(PCS_C_dur/h));
%---

for trial=1:cycle_num
    
    %---Pick a cycle length
    cycle_length=round(CS_dur(trial)/h);
    %---
    
    %--initialize values for timer and CS
    P=zeros(1,round(CS_dur(trial)/h));
    x=zeros(1,round(CS_dur(trial)/h));
    %--
    
    %---Conditions to retrieve rate and V from memory
    if CS_dur(trial)==CS_L_dur
        A=CS(2).rate(counterCS_L);
        V=CS(2).V(counterCS_L);
        counterCS_L=counterCS_L+1; % update counter for long CS
    elseif CS_dur(trial)==PCS_L_dur
        A=CS(2).rate(counterCS_L);
        V=CS(2).V(counterCS_L);
    elseif CS_dur(trial)==CS_S_dur
        A=CS(1).rate(counterCS_S);
        V=CS(1).V(counterCS_S);
        counterCS_S=counterCS_S+1; % update counter for short CS
    elseif CS_dur(trial)==PCS_S_dur
        A=CS(1).rate(counterCS_S);
        V=CS(1).V(counterCS_S);
    else
        A=(CS(1).rate(counterCS_S)+CS(2).rate(counterCS_L))/2;
        V=(CS(1).V(counterCS_S)+CS(2).V(counterCS_L))/2;
    end
        
    
    for t=1:cycle_length
        
        counterDDM=counterDDM+1; % update counter for random process in DDM
        
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
        CR(trial,t)=x(t)*V;
        %---
        
    end
    
    %---Slope Correction
    A=A+A*alpha_t*(1-P(t))/P(t); % realistic correction rule, never fully converges. Only updates in rewarded trials.
    %---
    
    %---V update
    V=RW(V,alpha_E,x(t),H,A,P(t));
    %---
    
    if CS_dur(trial)==CS_S_dur
        CS(1).rate(counterCS_S)=A;
        CS(1).V(counterCS_S)=V;
    elseif CS_dur(trial)==CS_L_dur 
        CS(2).rate(counterCS_L)=A;
        CS(2).V(counterCS_L)=V;
    end
    
end

CR_S_peak_mean=CR(CS_dur==PCS_S_dur,:);
CR_S_peak_mean=mean(CR_S_peak_mean(60:end,:),1);
CR_S_peak_mean=CR_S_peak_mean(1:round(50/h));

CR_L_peak_mean=CR(CS_dur==PCS_L_dur,:);
CR_L_peak_mean=mean(CR_L_peak_mean(60:end,:),1);
CR_L_peak_mean=CR_L_peak_mean(1:round(50/h));

CR_C_peak_mean=CR(CS_dur==PCS_C_dur,:);
CR_C_peak_mean=mean(CR_C_peak_mean(60:end,:),1);
CR_C_peak_mean=CR_C_peak_mean(1:round(50/h));

figure
plot(h:h:length(CR_S_peak_mean)*h,CR_S_peak_mean,h:h:length(CR_L_peak_mean)*h,CR_L_peak_mean,h:h:length(CR_C_peak_mean)*h,CR_C_peak_mean,'LineWidth',6)
PlotProperties
xlabel('time (s)')
ylabel('response strength')
legend('short','long','compound')

figure
[M_S,I_S]=max(CR_S_peak_mean);
[M_L,I_L]=max(CR_L_peak_mean);
[M_C,I_C]=max(CR_C_peak_mean);
plot((1:length(CR_S_peak_mean))/I_S,CR_S_peak_mean/M_S,(1:length(CR_L_peak_mean))/I_L,CR_L_peak_mean/M_L,(1:length(CR_C_peak_mean))/I_C,CR_C_peak_mean/M_C,'LineWidth',6)
xlabel('proportion of peak time')
ylabel('norm. response strength')
legend('short','long','compound')
PlotProperties

save('TempAvgSim','M_S','I_S','M_L','I_L','M_C','I_C','CR_S_peak_mean','CR_L_peak_mean','CR_C_peak_mean','h')
clear