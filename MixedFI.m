% Mixed FI

clear
close all

%---Parameters
h=0.01;
tau_x=50;
alpha_t=0.2;
alpha_E=0.1;
mu=1;
sigma=0.425;
m=0.2;
H=30;
%---

CS_dur=75; % (s)
US_A=15; % (s)
cycle_num=400;

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
CS_A=ones(1,round(CS_dur/h));
CS_B=ones(1,round(CS_dur/h));
%---

%---Initialize counter for DDM noise
counterDDM=0;
%---

%---CR
CR=NaN(cycle_num,round(CS_dur/h));
%---

%---Trial short or long
S_L=NaN(1,cycle_num);
S_L(1:round(length(S_L)/2))=round(CS_dur/h);
S_L(round(length(S_L)/2)+1:end)=round(US_A/h);
S_L=S_L(randperm(cycle_num));
%---

for trial=1:cycle_num
    
    %---Pick a cycle length
    cycle_length=S_L(trial);
    %---
    
    %--initialize values for timer and CS
    P_A=zeros(1,round(CS_dur/h));
    P_B=zeros(1,round(CS_dur/h));
    x_A=zeros(1,round(CS_dur/h));
    x_B=zeros(1,round(CS_dur/h));
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
        
        %---Timer selection
        P_att=min((P_A(t+1)-1)^2,(P_B(t+1)-1)^2);
        US_A_att=1*((P_A(t+1)-1)^2==P_att);
        US_B_att=1*((P_B(t+1)-1)^2==P_att);
        %---
        
        %---CR
        CR(trial,t)=(x_A(t)*V_A(trial)*US_A_att+x_B(t)*V_B(trial)*US_B_att);
        %---
        
    end
    
    %---Slope Correction
    A_A=A_A+US_A_att*A_A*alpha_t*(1-P_A(t))/P_A(t); % realistic correction rule, never fully converges. Only updates in rewarded trials.
    A_B=A_B+US_B_att*A_B*alpha_t*(1-P_B(t))/P_B(t); % realistic correction rule, never fully converges. Only updates in rewarded trials.
    %---
    
    %---V update
    V_A(trial+1)=RW(V_A(trial),alpha_E,US_A_att*x_A(t),H,A_A,P_A(t));
    V_B(trial+1)=RW(V_B(trial),alpha_E,US_B_att*x_B(t),H,A_B,P_B(t));
    %---
    
end

CR_mean=nanmean(CR(100:end,:),1);
plot(h:h:length(CR_mean)*h,CR_mean,'LineWidth',6)
title('model')
xlabel('time (s)')
ylabel('response strength')
axis([0 80 0 2])
PlotProperties

save('MixedFISim','CR_mean','h')
clear