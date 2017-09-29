% Variable Interval vs Fixed Interval peak procedure

clear
close all

%---Parameters
h=0.01;
tau_x=50;
alpha_t=0.1;
alpha_E=0.1;
mu=1;
sigma=0.3;
m=0.2;
H=40;
%---

%% VI

%--Create a random vector with trial durations
CS_dur=15:1:45; % (ms)
CS_dur=repmat(CS_dur,1,60);
CS_dur=[CS_dur 3*45*ones(1,round(0.25*length(CS_dur)))];
CS_dur=CS_dur(randperm(length(CS_dur)));
%---

cycle_num=length(CS_dur);

%---DDM constants
A_A=1*10^(-3);
N_A=normrnd(0,1,ceil(3*45/h)*cycle_num,1); % noise for CS A
%---

%---Associative strength and alpha
V_A=zeros(1,cycle_num); % initialize vector to store V values for time steps
%---

%---Initialize counter for DDM noise
counterDDM=0;
%---

%---CR
CR_VI=NaN(cycle_num,round(3*45/h));
%---

%---Trial short or long

%---

for trial=1:cycle_num
    
    %---Pick a cycle length
    cycle_length=round(CS_dur(trial)/h);
    %---
    
    %---Setting the presence or absence of CS and US
    CS_A=ones(1,cycle_length);
    %---
    
    %--initialize values for timer and CS
    P_A=zeros(1,cycle_length);
    x_A=zeros(1,cycle_length);
    %--
    
    for t=1:cycle_length
        
        counterDDM=counterDDM+1; % update counter for random process in DDM
        
        % min will take the minimum value: either DDM result or 3. This
        % caps the value of integrator at 3.
        P_A(t+1)=CS_A(t)*min(DDM( P_A(t), A_A, h, m, N_A(counterDDM) ), 3);
        
        % max ensures the minimum value the accumulator can reach is
        % 1*10^-6. This avoids division by zero later.
        P_A(t+1)=max(P_A(t+1), 1*10^(-3));
        
        %---Element (RBFs)
        x_A(t)=CStrace(P_A(t+1),mu,sigma,tau_x,CS_A(t),x_A(t),h);
        %---
        
        %---CR
        CR_VI(trial,t)=(x_A(t)*V_A(trial));
        %---
        
    end
    
    %---Slope Correction
    A_A=A_A+A_A*alpha_t*(1-P_A(t))/P_A(t); % realistic correction rule, never fully converges. Only updates in rewarded trials.
    %---
    
    %---V update
    V_A(trial+1)=RW(V_A(trial),alpha_E,x_A(t),H,A_A,P_A(t));
    %---
    
end

CR_VI_peak_mean=CR_VI(CS_dur==3*45,:);
CR_VI_peak_mean=mean(CR_VI_peak_mean(61:end,:),1);
CR_VI_peak_mean=CR_VI_peak_mean(1:round(80/h));

%% FI

%--Create a random vector with trial durations
CS_dur=30; % (ms)
CS_dur=repmat(CS_dur,1,500);
CS_dur=[CS_dur 3*30*ones(1,round(0.25*length(CS_dur)))];
CS_dur=CS_dur(randperm(length(CS_dur)));
%---

cycle_num=length(CS_dur);

%---DDM constants
A_A=1*10^(-3);
N_A=normrnd(0,1,ceil(3*30/h)*cycle_num,1); % noise for CS A
%---

%---Associative strength and alpha
V_A=zeros(1,cycle_num); % initialize vector to store V values for time steps
%---

%---Initialize counter for DDM noise
counterDDM=0;
%---

%---CR
CR_FI=NaN(cycle_num,round(3*30/h));
%---

%---Trial short or long

%---

for trial=1:cycle_num
    
    %---Pick a cycle length
    cycle_length=round(CS_dur(trial)/h);
    %---
    
    %---Setting the presence or absence of CS and US
    CS_A=ones(1,cycle_length);
    %---
    
    %--initialize values for timer and CS
    P_A=zeros(1,cycle_length);
    x_A=zeros(1,cycle_length);
    %--
    
    for t=1:cycle_length
        
        counterDDM=counterDDM+1; % update counter for random process in DDM
        
        % min will take the minimum value: either DDM result or 3. This
        % caps the value of integrator at 3.
        P_A(t+1)=CS_A(t)*min(DDM( P_A(t), A_A, h, m, N_A(counterDDM) ), 3);
        
        % max ensures the minimum value the accumulator can reach is
        % 1*10^-6. This avoids division by zero later.
        P_A(t+1)=max(P_A(t+1), 1*10^(-3));
        
        %---Element (RBFs)
        x_A(t)=CStrace(P_A(t+1),mu,sigma,tau_x,CS_A(t),x_A(t),h);
        %---
        
        %---CR
        CR_FI(trial,t)=(x_A(t)*V_A(trial));
        %---
        
    end
    
    %---Slope Correction
    A_A=A_A+A_A*alpha_t*(1-P_A(t))/P_A(t); % realistic correction rule, never fully converges. Only updates in rewarded trials.
    %---
    
    %---V update
    V_A(trial+1)=RW(V_A(trial),alpha_E,x_A(t),H,A_A,P_A(t));
    %---
    
end

CR_FI_peak_mean=CR_FI(CS_dur==3*30,:);
CR_FI_peak_mean=mean(CR_FI_peak_mean(61:end,:),1);
CR_FI_peak_mean=CR_FI_peak_mean(1:round(80/h));

figure
plot(h:h:length(CR_VI_peak_mean)*h,CR_VI_peak_mean,h:h:length(CR_FI_peak_mean)*h,CR_FI_peak_mean,'LineWidth',6)
legend({'VI','FI'},'Box','off')
title('model')
xlabel('time (s)')
ylabel('response strength')
PlotProperties

% scalar plot
[M_FI,I_FI]=max(CR_FI_peak_mean);
[M_VI,I_VI]=max(CR_VI_peak_mean);
figure
plot((1:length(CR_VI_peak_mean))/I_VI,CR_VI_peak_mean/M_VI,(1:length(CR_FI_peak_mean))/I_FI,CR_FI_peak_mean/M_FI,'LineWidth',6)
legend({'VI','FI'},'Box','off')
title('model')
xlabel('proportion of peak time')
ylabel('response strength')
PlotProperties

save('VI_FI_Sim','M_FI','I_FI','M_VI','I_VI','CR_VI_peak_mean','CR_FI_peak_mean','h')
clear
