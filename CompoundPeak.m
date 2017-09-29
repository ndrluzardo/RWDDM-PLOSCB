% Compound peak procedure

clear
close all
%% 
% Two CSs are used: CS A (conditioned alone) and CS B (conditioned alone). 
% Then both CSs are presented as a compound in a trial that lasts 3 times
% longer.
% 
% This model will select the timer nearer to the threshold to pay attention. 
%% Trial Structure
% The cycle is composed of CS, US and ITI in this order as follows:
% 
%  
% 
% Choose the duration of each CS (in milliseconds) using the variables below:

CS_dur=50; % (ms)
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
%% Phase 2: Compound conditioning
% The two CSs will be presented in compound. No reinforcement will be given 
% and there will be no update to the V values or timer slopes. Normal
% conditioning trials will be insterspersed.
% 
% Enter the number of trials in the variable _cycle_num_ below.

cycle_num=700;

percent_peak=0.25;
percent_compound=1/3; % percentage of compound trials

peak_trials=randperm(cycle_num); % actual trial number for each peak trial
peak_trials=peak_trials(1:round(cycle_num*percent_peak));

limits=linspace(1,length(peak_trials),4); % indices that divide peak_trials vector in three

compound_trials=peak_trials(1:round(limits(2))); % actual trial number for each compound trial
single_trials_A=peak_trials(round(limits(2))+1:round(limits(3))); % actual trial number for each single peak trial
single_trials_B=peak_trials(round(limits(3))+1:end);
%% 
% 

cycle_length=CS_dur/h;

%---DDM constants
N_A=normrnd(0,1,ceil(CS_dur*3/h)*cycle_num,1); % noise for CS A
N_B=normrnd(0,1,ceil(CS_dur*3/h)*cycle_num,1); % noise for CS B
%---

%---Initialize counters for DDM noise and CRs
counterDDM=0;
counterCRAB=0;
counterCRAorB=0;
%---

%---Carry over V values from previous phase
V_A=V_A(end);
V_B=V_B(end);
%---

%---Initialize CR vector
CR_AB=zeros(length(compound_trials),cycle_length*3);
CR_AorB=zeros(length(single_trials_A)+length(single_trials_B),cycle_length*3);
%---

for trial=1:cycle_num
    
    if ismember(trial,compound_trials)==1 % compound peak trial
        counterCRAB=counterCRAB+1;
        
        %---Setting the presence or absence of CS and US
        CS_A=ones(1,cycle_length*3);
        CS_B=ones(1,cycle_length*3);
        %---
        
        %--initialize values for timer and CS
        P_A=zeros(1,cycle_length*3);
        P_B=zeros(1,cycle_length*3);
        x_A=zeros(1,cycle_length*3);
        x_B=zeros(1,cycle_length*3);
        %--
        
        for t=1:cycle_length*3
            
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
        %Att=[P_A(t+1)>P_B(t+1) P_B(t+1)>P_A(t+1)]*1;
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
    elseif ismember(trial,single_trials_A)==1 % peak trial CSA
        counterCRAorB=counterCRAorB+1;
        
        %---Setting the presence or absence of CS and US
        CS_A=ones(1,cycle_length*3);
        %---
        
        %--initialize values for timer and CS
        P_A=zeros(1,cycle_length*3);
        x_A=zeros(1,cycle_length*3);
        %--
        
        for t=1:cycle_length*3
            
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
            CR_AorB(counterCRAorB,t)=x_A(t)*V_A;
            %---
            
        end
        
    elseif ismember(trial,single_trials_B)==1 % peak trial CSB
        counterCRAorB=counterCRAorB+1;
        
        %---Setting the presence or absence of CS and US
        CS_B=ones(1,cycle_length*3);
        %---
        
        %--initialize values for timer and CS
        P_B=zeros(1,cycle_length*3);
        x_B=zeros(1,cycle_length*3);
        %--
        
        for t=1:cycle_length*3
            
        counterDDM=counterDDM+1; % update counter for random process in DDM
        
        % min will take the minimum value: either DDM result or 3. This
        % caps the value of integrator at 3.
        P_B(t+1)=CS_B(t)*min(DDM( P_B(t), A_B, h, m, N_B(counterDDM) ), 3);
        
        % max ensures the minimum value the accumulator can reach is
        % 1*10^-6. This avoids division by zero later.
        P_B(t+1)=max(P_B(t+1), 1*10^(-6));
        
        
        %---Element (RBFs)
        x_B(t)=CStrace(P_B(t+1),mu,sigma,tau_x,CS_B(t),x_B(t),h);
        %---
        
        %---CR
        CR_AorB(counterCRAorB,t)=(x_B(t)*V_B);
        %---
        
        end
    else % reinforced trial with single CSA or CSB
        
        %---Setting the presence or absence of CS and US
        CS_A=ones(1,cycle_length);
        CS_B=ones(1,cycle_length);
        %---
        
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
            P_A(t+1)=max(P_A(t+1), 1*10^(-6));
            P_B(t+1)=max(P_B(t+1), 1*10^(-6));
            
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
        V_A=RW(V_A,alpha_E,x_A(t),H,A_A,P_A(t));
        V_B=RW(V_B,alpha_E,x_B(t),H,A_B,P_B(t));
        %---
    end
end

figure
CR_AorB_mean=mean(CR_AorB,1);
CR_AB_mean=mean(CR_AB,1);
plot(h:h:length(CR_AorB_mean)*h,CR_AorB_mean,h:h:length(CR_AB_mean)*h,CR_AB_mean,'LineWidth',6)
xlabel('time (s)')
ylabel('response strength')
title('model')
legend({'CR A and B','CR AB'},'Box','off')
axis([0 150 0 0.15])
PlotProperties

figure
CR_AorB_perc=CR_AorB_mean./max(CR_AorB_mean);
CR_AB_perc=CR_AB_mean./max(CR_AB_mean);
plot(h:h:length(CR_AorB_perc)*h,CR_AorB_perc,h:h:length(CR_AB_perc)*h,CR_AB_perc,'LineWidth',6)
xlabel('time (s)')
ylabel('norm. response strength')
title('model')
legend({'CR A and B','CR AB'},'Box','off')
PlotProperties

save('CompoundPeakSim','CR_AorB_mean','CR_AB_mean','CR_AorB_perc','CR_AB_perc','h')

PeakCSAorB=find(CR_AorB_mean==max(CR_AorB_mean))*h/1000;
PeakCSAB=find(CR_AB_mean==max(CR_AB_mean))*h/1000;