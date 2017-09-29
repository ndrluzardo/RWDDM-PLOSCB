% Acquistion and Extinction

clear
close all

CSA_dur=5; % (s)

%---Parameters for a LTDDM with RW-PH
h=0.01;
tau_x=50;
alpha_t=0.1;
alpha_E=0.1;
mu=1;
sigma=0.30;
m=0.15;
H=4;
%---

acq_num=80;

ext_num=100;

cycle_num=acq_num+ext_num;

%---DDM constants
N=normrnd(0,1,ceil(CSA_dur/h)*cycle_num,1); % noise for CS
%---

cycle_length=round(CSA_dur/h);

%---Associative strength and A
V=zeros(1,cycle_num); % initialize vector to store V values for time steps
A=zeros(1,cycle_num);
A(1)=1*10^(-3);
%---

%---Setting the presence or absence of CS and US
CS=ones(1,cycle_length); % CS is 1 only during CS presentation
acq_trial=zeros(1,cycle_num); % set up an acquisition vector
acq_trial(1,1:acq_num)=1; % acquisition vector is 1 if trial is acquisition, 0 otherwise
%---

%---Initialize counter for DDM noise
counterDDM=0;
%---


for trial=1:cycle_num
    
    %--initialize values for timer and CS
    P=zeros(1,cycle_length);
    x=zeros(1,cycle_length);
    %--
    
    for t=1:cycle_length
        
        counterDDM=counterDDM+1; % update counter for random process in DDM
        
        % min will take the minimum value: either DDM result or 3. This
        % caps the value of integrator at 3.
        P(t+1)=CS(t)*min(DDM( P(t), A(trial), h, m, N(counterDDM) ), 3);
        
        % max ensures the minimum value the accumulator can reach is
        % 1*10^-6. This avoids division by zero later.
        P(t+1)=max(P(t+1), 1*10^(-3));
        
        %---Element (RBFs)
        x(t+1)=CStrace(P(t+1),mu,sigma,tau_x,CS(t),x(t),h);
        %---
    end
    
    %---V update
    V(trial+1)=RW( V(trial),alpha_E,x(t),acq_trial(trial)*H,A(trial),P(t));
    %---
    
    %---Slope Correction
    A(trial+1)=A(trial)+A(trial)*alpha_t*(1-P(t))/P(t); % realistic correction rule, never fully converges. Only updates in rewarded trials.
    %---
    
end

%% Reacquisition

cycle_num=acq_num;

%---Associative strength and A
Vreac=zeros(1,cycle_num); % initialize vector to store V values for time steps
Areac=zeros(1,cycle_num);
Areac(1)=A(end);
%---

%---DDM constants
N=normrnd(0,1,ceil(CSA_dur/h)*cycle_num,1); % noise for CS
%---
%---Initialize counter for DDM noise
counterDDM=0;
%---

for trial=1:cycle_num
    
    %--initialize values for timer and CS
    P=zeros(1,cycle_length);
    x=zeros(1,cycle_length);
    %--
    
    for t=1:cycle_length
        
        counterDDM=counterDDM+1; % update counter for random process in DDM
        
        % min will take the minimum value: either DDM result or 3. This
        % caps the value of integrator at 3.
        P(t+1)=CS(t)*min(DDM( P(t), Areac(trial), h, m, N(counterDDM) ), 3);
        
        % max ensures the minimum value the accumulator can reach is
        % 1*10^-6. This avoids division by zero later.
        P(t+1)=max(P(t+1), 1*10^(-3));
        
        %---Element (RBFs)
        x(t+1)=CStrace(P(t+1),mu,sigma,tau_x,CS(t),x(t),h);
        %---
    end
    
    %---V update
    Vreac(trial+1)=RW( Vreac(trial),alpha_E,x(t),acq_trial(trial)*H,Areac(trial),P(t));
    %---
    
    %---Slope Correction
    Areac(trial+1)=Areac(trial)+Areac(trial)*alpha_t*(1-P(t))/P(t); % realistic correction rule, never fully converges. Only updates in rewarded trials.
    %---
    
end

%% figures


% The next figure shows activation of the CS A representation:
% figure
% plot(h/1000:h/1000:length(x)*h/1000,x,'LineWidth',6)
% PlotProperties
% xlabel('time (sec)')
% ylabel('CS A activation')
% 
% This plot shows the evolution of the V value over trials:
figure
plot(1:length(V),V,'LineWidth',6)
PlotProperties
ylim([0,1])
title('model')
xlabel('trials')
ylabel('associative strength')

% plots the evolution of A over trials
figure
plot(1:length(A),A,'LineWidth',6)
PlotProperties
%yticks(0:
title('model')
xlabel('trials')
ylabel('slope A')

% plots acquisition and reacquisition
figure
plot(1:acq_num,V(1:acq_num),1:acq_num,Vreac(1:acq_num),'LineWidth',6)
PlotProperties
title('model')
xlabel('trials')
ylabel('associative strength')
legend({'acquisition','reacquisition'},'Location','southeast','Box','off')

save('AcqReacqSimData','V','A','acq_num','Vreac')
clear
