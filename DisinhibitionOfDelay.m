%% Disinhibition of Delay

clear
close all
clc
%% 
% Two CSs are used: CS A (conditioned alone) and CS B (conditioned alone). 
% Then both CSs are presented as a compound.
% 
% This model will select the timer nearer to the threshold to pay attention. 
%% Trial Structure
% The cycle is composed of CS, US and ITI in this order as follows:
% 
%  
% 
% Choose the duration of each CS (in milliseconds) using the variables below:

CS_dur=30; % (ms)

%% Parameters

h=0.01;
tau_x=50;
alpha_t=0.75;
alpha_E=0.1;
mu=1;
sigma=0.18;
m=0.25;
H=5;

%% Phase 1: Conditioning CS A and CS B
% CS A and CS B will be separately conditioned (A+, B+).
% 
% Enter the number of conditioning trials in the variable _cycle_num_ below.

cycle_num=100;
%% 

cycle_length=round(CS_dur/h);
CSA_dur=CS_dur;
CSB_dur=CS_dur;

%---DDM constants
A_A=1*10^(-3);
A_B=1*10^(-3);
N_A=normrnd(0,1,ceil(CS_dur/h)*cycle_num,1); % noise for CS A
N_B=normrnd(0,1,ceil(CS_dur/h)*cycle_num,1); % noise for CS B
%---

%---Associative strength and alpha
V_A=zeros(1,cycle_num); % initialize vector to store V values for time steps
V_B=zeros(1,cycle_num); % initialize vector to store V values for time steps
%---

%---Setting the presence or absence of CS and US
CS_A=ones(1,cycle_length);
CS_B=ones(1,cycle_length);
%---

%---Initialize counter for DDM noise
counterDDM=0;
%---
for trial=1:cycle_num
    
    %--initialize values for timer and CS
    P_A=zeros(1,cycle_length);
    P_B=zeros(1,cycle_length);
    x_A=zeros(1,cycle_length);
    x_B=zeros(1,cycle_length);
    %--
    
    for t=1:cycle_length
        
        counterDDM=counterDDM+1; % update counter for random process in DDM
        
        % min will take the minimum value: either DDM result or 3. This
        % caps the value of integrator at 3.
        P_A(t+1)=CS_A(t)*min(DDM( P_A(t), A_A, h, m, N_A(counterDDM) ), 3);
        P_B(t+1)=CS_B(t)*min(DDM( P_B(t), A_B, h, m, N_B(counterDDM) ), 3);
        
        % max ensures the minimum value the accumulator can reach is
        % 1*10^-6. This avoids division by zero later.
        P_A(t+1)=max(P_A(t+1), 1*10^(-3));
        P_B(t+1)=max(P_B(t+1), 1*10^(-3));
        
        %---Element (RBFs)
        x_A(t)=CStrace(P_A(t+1),mu,sigma,tau_x,CS_A(t),x_A(t),h);
        x_B(t)=CStrace(P_B(t+1),mu,sigma,tau_x,CS_B(t),x_B(t),h);
        %---
        
    end
    
    %---Slope Correction
    A_A=A_A+A_A*alpha_t*(1-P_A(t))/P_A(t); % realistic correction rule, never fully converges. Only updates in rewarded trials.
    A_B=A_B+A_B*alpha_t*(1-P_B(t))/P_B(t); % realistic correction rule, never fully converges. Only updates in rewarded trials.
    %---
    
    %---V update
    V_A(trial+1)=RW(V_A(trial),alpha_E,x_A(t),H,A_A,P_A(t));
    V_B(trial+1)=RW(V_B(trial),alpha_E,x_B(t),H,A_B,P_B(t));
    %---
    
end
%% 
% Here we plot the evolution of V values over trials for both CSs

% plot(1:length(V_A),V_A,1:length(V_B),V_B,'LineWidth',3)
% PlotProperties
% xlabel('trials')
% ylabel('associative strength')
% legend('CS A','CS B')
%% 
% Plots CS representations:

% plot(1:length(x_A),x_A,1:length(x_B),x_B,'LineWidth',3)
% PlotProperties
% legend('CS A','CS B')
% xlabel('time')
% ylabel('CS activation')
%% Phase 2: Compound conditioning
% The two CSs will be presented in compound. No reinforcement will be given 
% and there will be no update to the V values or timer slopes. Normal
% conditioning trials will be insterspersed.
% 
% Enter the number of trials in the variable _cycle_num_ below.

cycle_num=500;

percent_compound=0.5; % percentage of compound trials

compound_trials=randperm(cycle_num); % actual trial number for each compound trial
compound_trials=compound_trials(1:round(cycle_num*percent_compound));
%% 
% 

cycle_length=CS_dur/h;

%---DDM constants
N_A=normrnd(0,1,ceil(CS_dur/h)*cycle_num,1); % noise for CS A
N_B=normrnd(0,1,ceil(CS_dur/h)*cycle_num,1); % noise for CS B
%---

%---Initialize counters for DDM noise and CRs
counterDDM=0;
counterCRAB=0;
counterCRAorB=0;
%---

%---Initialize CR vector
CR_AB=zeros(length(compound_trials),cycle_length);
CR_A=zeros(cycle_num-length(compound_trials),cycle_length);
CR_B=zeros(cycle_num-length(compound_trials),cycle_length);
%---

%---Carry over V values from previous phase
V_A=V_A(end);
V_B=V_B(end);
%---

for trial=1:cycle_num
    
    %--initialize values for timer and CS
    P_A=zeros(1,cycle_length);
    P_B=zeros(1,cycle_length);
    x_A=zeros(1,cycle_length);
    x_B=zeros(1,cycle_length);
    %--
    
    if ismember(trial,compound_trials)==1
        counterCRAB=counterCRAB+1;
        
        for t=1:cycle_length
        
        counterDDM=counterDDM+1; % update counter for random process in DDM
        
        % min will take the minimum value: either DDM result or 3. This
        % caps the value of integrator at 3.
        P_A(t+1)=CS_A(t)*min(DDM( P_A(t), A_A, h, m, N_A(counterDDM) ), 3);
        P_B(t+1)=CS_B(t)*min(DDM( P_B(t), A_B, h, m, N_B(counterDDM) ), 3);
        
        % max ensures the minimum value the accumulator can reach is
        % 1*10^-6. This avoids division by zero later.
        P_A(t+1)=max(P_A(t+1), 1*10^(-3));
        P_B(t+1)=max(P_B(t+1), 1*10^(-3));
        
        %--- attention to timer
        Att=[A_A>A_B A_B>A_A]*1;
        %---
        
        %---Element (RBFs)
        x_A(t)=CStrace(P_A(t+1),mu,sigma,tau_x,CS_A(t),x_A(t),h);
        x_B(t)=CStrace(P_B(t+1),mu,sigma,tau_x,CS_B(t),x_B(t),h);
        %---
        
        %---CR
        CR_AB(counterCRAB,t)=(x_A(t)*V_A*Att(1)+x_B(t)*V_B*Att(2));
        %---
        
        end
    else
        counterCRAorB=counterCRAorB+1;
        
        for t=1:cycle_length
            
            counterDDM=counterDDM+1; % update counter for random process in DDM
            
            % min will take the minimum value: either DDM result or 3. This
            % caps the value of integrator at 3.
            P_A(t+1)=CS_A(t)*min(DDM( P_A(t), A_A, h, m, N_A(counterDDM) ), 3);
            P_B(t+1)=CS_B(t)*min(DDM( P_B(t), A_B, h, m, N_B(counterDDM) ), 3);
            
            % max ensures the minimum value the accumulator can reach is
            % 1*10^-6. This avoids division by zero later.
            P_A(t+1)=max(P_A(t+1), 1*10^(-3));
            P_B(t+1)=max(P_B(t+1), 1*10^(-3));
            
            %---Element (RBFs)
            x_A(t)=CStrace(P_A(t+1),mu,sigma,tau_x,CS_A(t),x_A(t),h);
            x_B(t)=CStrace(P_B(t+1),mu,sigma,tau_x,CS_B(t),x_B(t),h);
            %---
            
            %---CR
            CR_A(counterCRAorB,t)=x_A(t)*V_A(end);
            CR_B(counterCRAorB,t)=x_B(t)*V_B(end);
            %---
            
        end
        
        %---Slope Correction
        A_A=A_A+A_A*alpha_t*(1-P_A(t))/P_A(t); % realistic correction rule, never fully converges. Only updates in rewarded trials.
        A_B=A_B+A_B*alpha_t*(1-P_B(t))/P_B(t); % realistic correction rule, never fully converges. Only updates in rewarded trials.
        %---
        
        %---V update
        V_A=RW(V_A,alpha_E,x_A(t),H,A_A,P_A(t));
        V_B=RW(V_B,alpha_E,x_B(t),H,A_B,P_B(t));
        %---
        
        
    end
end

CR_AandB_mean=mean(vertcat(CR_A,CR_B),1);
CR_AB_mean=mean(CR_AB,1);
plot(h:h:length(CR_AandB_mean)*h,CR_AandB_mean,h:h:length(CR_AB_mean)*h,CR_AB_mean,'LineWidth',6)
title('model')
xlabel('time (s)')
ylabel('response strength')
legend({'CR A and B','CR AB'},'Location','northwest','Box','off')
PlotProperties

save('DisinhDelaySim','CR_AandB_mean','CR_AB_mean','h')