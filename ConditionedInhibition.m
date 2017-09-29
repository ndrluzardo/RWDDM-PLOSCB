% Conditioned Inhibition

clear
close all

%% phase 1 - Acquisition only with short and long excitors and inhibitors

%---cycles
cycle_num=1200;
%---

%---One CS, two durations
CS_dur_L=30;
CS_dur_S=10;
%---

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
CS(1).name='Excitor Short';
CS(2).name='Inhibitor Short';
CS(3).name='Excitor Long';
CS(4).name='Inhibitor Long';
%-A
CS(1).A.ph1=zeros(cycle_num,1);
CS(1).A.ph1(1)=1*10^(-3);

CS(2).A.ph1=zeros(cycle_num,1);
CS(2).A.ph1(1)=1*10^(-3);

CS(3).A.ph1=zeros(cycle_num,1);
CS(3).A.ph1(1)=1*10^(-3);

CS(4).A.ph1=zeros(cycle_num,1);
CS(4).A.ph1(1)=1*10^(-3);
%-V
CS(1).V.ph1=zeros(cycle_num,1);
CS(2).V.ph1=zeros(cycle_num,1);
CS(3).V.ph1=zeros(cycle_num,1);
CS(4).V.ph1=zeros(cycle_num,1);
%---

%---Trial short or long
S_L=NaN(3,cycle_num);
S_L(1,1:round(length(S_L)/2))=round(CS_dur_L/h);
S_L(2,1:round(length(S_L)/2))=3;
S_L(3,1:round(length(S_L)/2))=4;
S_L(1,round(length(S_L)/2)+1:end)=round(CS_dur_S/h);
S_L(2,round(length(S_L)/2)+1:end)=1;
S_L(3,round(length(S_L)/2)+1:end)=2;
S_L=S_L(:,randperm(cycle_num));
%---

%---Trial excitation or inhibition
E_I=NaN(1,cycle_num);
E_I(1:round(length(E_I)/2))=1; % excitation is 1
E_I(round(length(E_I)/2)+1:end)=0; % inhibition is 0
E_I=E_I(randperm(cycle_num));
%---

%---DDM constants
N1=normrnd(0,1,ceil(CS_dur_L/h*cycle_num),1); % noise for CS A
N2=normrnd(0,1,ceil(CS_dur_L/h*cycle_num),1); % noise for CS A
%---

counterDDM=0;

for trial=1:cycle_num
    
    cycle_length=S_L(1,trial);
    
    %--initialize values for timer and CS. 1=excitor, 2=inhibitor
    P1=zeros(1,cycle_length);
    P2=zeros(1,cycle_length);

    x1=zeros(1,cycle_length);
    x2=zeros(1,cycle_length);
    %--
    
    for t=1:cycle_length
        
        counterDDM=counterDDM+1; % update counter for random process in DDM
        
        A1=CS(S_L(2,trial)).A.ph1(trial); % excitor
        A2=CS(S_L(3,trial)).A.ph1(trial); % inhibitor
        V1=CS(S_L(2,trial)).V.ph1(trial); % excitor
        V2=CS(S_L(3,trial)).V.ph1(trial); % inhibitor
        
        % min will take the minimum value: either DDM result or 3. This
        % caps the value of integrator at 3.
        P1(t+1)=min(DDM( P1(t), A1, h, m, N1(counterDDM) ), 3);
        P2(t+1)=min(DDM( P2(t), A2, h, m, N2(counterDDM) ), 3);

        % max ensures the minimum value the accumulator can reach is
        % 1*10^-6. This avoids division by zero later.
        P1(t+1)=max(P1(t+1), 1*10^(-3));
        P2(t+1)=max(P2(t+1), 1*10^(-3));
        
        %---Element (RBFs)
        x1(t)=CStrace(P1(t+1),mu,sigma,tau_x,1,x1(t),h);
        x2(t)=(1-E_I(trial))*CStrace(P2(t+1),mu,sigma,tau_x,1,x2(t),h);
        %---
        
    end
    
    %---Slope Correction
    A1=A1+A1*alpha_t*(1-P1(t+1))/P1(t+1); % realistic correction rule, never fully converges. Only updates in rewarded trials.
    A2=A2+(1-E_I(trial))*A2*alpha_t*(1-P2(t+1))/P2(t+1);
    %---
    
    %---V update
    V1=RW([V1 V2],alpha_E,[x1(t) x2(t)],E_I(trial)*H,A1,P1(t+1));
    V2=RW([V2 V1],alpha_E,(1-E_I(trial))*[x2(t) x1(t)],E_I(trial)*H,A2,P2(t+1));
    %---
    
    %---Move every A and V value 1 trial forward
    CS(1).A.ph1(trial+1)=CS(1).A.ph1(trial);
    CS(2).A.ph1(trial+1)=CS(2).A.ph1(trial);
    CS(3).A.ph1(trial+1)=CS(3).A.ph1(trial);
    CS(4).A.ph1(trial+1)=CS(4).A.ph1(trial);
    CS(1).V.ph1(trial+1)=CS(1).V.ph1(trial);
    CS(2).V.ph1(trial+1)=CS(2).V.ph1(trial);
    CS(3).V.ph1(trial+1)=CS(3).V.ph1(trial);
    CS(4).V.ph1(trial+1)=CS(4).V.ph1(trial);
    %---
    
    %---Now update only the As and Vs used in previous trial
    CS(S_L(2,trial)).A.ph1(trial+1)=A1; % excitor
    CS(S_L(3,trial)).A.ph1(trial+1)=A2; % inhibitor
    CS(S_L(2,trial)).V.ph1(trial+1)=V1; % excitor
    CS(S_L(3,trial)).V.ph1(trial+1)=V2; % inhibitor
    %---
    
end

%% Phase 2 - Mixed FI transfer excitor

%---cycles
cycle_num=600;
%---

%---structure arrays
CS(5).name='Transfer Excitor';
%-A
CS(5).A(1).ph1=zeros(cycle_num,1);
CS(5).A(1).ph1(1)=1*10^(-3);

CS(5).A(2).ph1=zeros(cycle_num,1);
CS(5).A(2).ph1(1)=1*10^(-3);

%-V
CS(5).V(1).ph1=zeros(cycle_num,1);
CS(5).V(2).ph1=zeros(cycle_num,1);
%---

%---Trial short or long: 0 for short, 1 for long
S_L=NaN(2,cycle_num);
S_L(1,1:round(length(S_L)/2))=round(CS_dur_L/h);
S_L(2,1:round(length(S_L)/2))=1;
S_L(1,round(length(S_L)/2)+1:end)=round(CS_dur_S/h);
S_L(2,round(length(S_L)/2)+1:end)=0;
S_L=S_L(:,randperm(cycle_num));
%---

%---DDM constants
N1=normrnd(0,1,ceil(CS_dur_L/h*cycle_num),1); % noise for CS A
N2=normrnd(0,1,ceil(CS_dur_L/h*cycle_num),1); % noise for CS A
%---


counterDDM=0;

for trial=1:cycle_num
    
    cycle_length=S_L(1,trial);
    
    %--initialize values for timer and CS. 1=excitor, 2=inhibitor
    P1=zeros(1,cycle_length);
    P2=zeros(1,cycle_length);

    x1=zeros(1,cycle_length);
    x2=zeros(1,cycle_length);
    %--
    
    for t=1:cycle_length
        
        counterDDM=counterDDM+1; % update counter for random process in DDM
        
        A1=CS(5).A(1).ph1(trial); % excitor
        A2=CS(5).A(2).ph1(trial); % inhibitor
        V1=CS(5).V(1).ph1(trial); % excitor
        V2=CS(5).V(2).ph1(trial); % inhibitor
        
        % min will take the minimum value: either DDM result or 3. This
        % caps the value of integrator at 3.
        P1(t+1)=min(DDM( P1(t), A1, h, m, N1(counterDDM) ), 3);
        P2(t+1)=min(DDM( P2(t), A2, h, m, N2(counterDDM) ), 3);

        % max ensures the minimum value the accumulator can reach is
        % 1*10^-6. This avoids division by zero later.
        P1(t+1)=max(P1(t+1), 1*10^(-3));
        P2(t+1)=max(P2(t+1), 1*10^(-3));
        
        %---Timer selection
        P_att=min((P1(t+1)-1)^2,(P2(t+1)-1)^2);
        US_A_att=1*((P1(t+1)-1)^2==P_att);
        US_B_att=1*((P2(t+1)-1)^2==P_att);
        %---
        
        %---Element (RBFs)
        x1(t)=CStrace(P1(t+1),mu,sigma,tau_x,1,x1(t),h);
        x2(t)=CStrace(P2(t+1),mu,sigma,tau_x,1,x2(t),h);
        %---
        
    end
    
    %---Slope Correction
    A1=A1+US_A_att*A1*alpha_t*(1-P1(t+1))/P1(t+1); % realistic correction rule, never fully converges. Only updates in rewarded trials.
    A2=A2+US_B_att*A2*alpha_t*(1-P2(t+1))/P2(t+1);
    %---
    
    %---V update
    V1=RW(V1,alpha_E,US_A_att*x1(t),H,A1,P1(t+1));
    V2=RW(V2,alpha_E,US_B_att*x2(t),H,A2,P2(t+1));
    %---
    
    %---Now update only the As and Vs used in previous trial
    CS(5).A(1).ph1(trial+1)=A1; % excitor
    CS(5).A(2).ph1(trial+1)=A2; % inhibitor
    CS(5).V(1).ph1(trial+1)=V1; % excitor
    CS(5).V(2).ph1(trial+1)=V2; % inhibitor
    %---
end

%% Test phase - peak procedure only trials

%---cycles
cycle_num=100;
%---

%---Trial short or long: 1 for short, 2 for long
S_L=NaN(1,cycle_num);
S_L(1,1:round(length(S_L)/2))=1;
S_L(1,round(length(S_L)/2)+1:end)=2;
S_L=S_L(:,randperm(cycle_num));
%---

%---DDM constants
N1=normrnd(0,1,ceil(3*CS_dur_L/h*cycle_num),1); % noise for CS A
N2=normrnd(0,1,ceil(3*CS_dur_L/h*cycle_num),1); % noise for CS A
N3=normrnd(0,1,ceil(3*CS_dur_L/h*cycle_num),1);
%---

counterDDM=0;

cycle_length=round(3*CS_dur_L/h);

%---Use an average for A and V, since they are not updated
A1=mean(CS(5).A(1).ph1(end-150:end)); % transfer excitor A1
A2=mean(CS(5).A(2).ph1(end-150:end)); % transfer excitor A2
V1=mean(CS(5).V(1).ph1(end-150:end)); % transfer excitor V1
V2=mean(CS(5).V(2).ph1(end-150:end)); % transfer excitor V2
A3(1)=mean(CS(2).A.ph1(end-150:end)); % inhibitor short A
V3(1)=mean(CS(2).V.ph1(end-150:end)); % inhibitor short V
A3(2)=mean(CS(4).A.ph1(end-150:end)); % inhibitor long A
V3(2)=mean(CS(4).V.ph1(end-150:end)); % inhibitor long V
%---

%---CR matrix
CR=NaN(cycle_num,cycle_length);
CR_exc=NaN(cycle_num,cycle_length);
%---

for trial=1:cycle_num
    
    %--initialize values for timer and CS. 1=excitor, 2=inhibitor
    P1=zeros(1,cycle_length);
    P2=zeros(1,cycle_length);
    P3=zeros(1,cycle_length);

    x1=zeros(1,cycle_length);
    x2=zeros(1,cycle_length);
    x3=zeros(1,cycle_length);
    %--
    
    
    for t=1:cycle_length
        
        counterDDM=counterDDM+1; % update counter for random process in DDM
        
        % min will take the minimum value: either DDM result or 3. This
        % caps the value of integrator at 3.
        P1(t+1)=min(DDM( P1(t), A1, h, m, N1(counterDDM) ), 3);
        P2(t+1)=min(DDM( P2(t), A2, h, m, N2(counterDDM) ), 3);
        P3(t+1)=min(DDM( P3(t), A3(S_L(trial)), h, m, N3(counterDDM) ), 3);
        
        % max ensures the minimum value the accumulator can reach is
        % 1*10^-6. This avoids division by zero later.
        P1(t+1)=max(P1(t+1), 1*10^(-3));
        P2(t+1)=max(P2(t+1), 1*10^(-3));
        P3(t+1)=max(P3(t+1), 1*10^(-3));
        
        %---Timer selection
        P_att=min((P1(t+1)-1)^2,(P2(t+1)-1)^2);
        US_A_att=1*((P1(t+1)-1)^2==P_att);
        US_B_att=1*((P2(t+1)-1)^2==P_att);
        %---
        
        %---Element (RBFs)
        x1(t)=CStrace(P1(t+1),mu,sigma,tau_x,1,x1(t),h);
        x2(t)=CStrace(P2(t+1),mu,sigma,tau_x,1,x2(t),h);
        x3(t)=CStrace(P3(t+1),mu,sigma,tau_x,1,x3(t),h);
        %---
        
        %---CR with inhibitor
        CR(trial,t)=(x1(t)*V1*US_A_att+x2(t)*V2*US_B_att+x3(t)*V3(S_L(trial)))...
            *((x1(t)*V1*US_A_att+x2(t)*V2*US_B_att+x3(t)*V3(S_L(trial)))>=0);
        %---
        
        %---CR without inhibitor
        CR_exc(trial,t)=(x1(t)*V1*US_A_att+x2(t)*V2*US_B_att)*((x1(t)*V1*US_A_att+x2(t)*V2*US_B_att)>=0);
        %---
        
    end
    
end

CR_mean_exc=mean(CR_exc(1:round(cycle_num/2),1:round(80/h)));
CR_mean_inh_short=mean(CR(S_L==1,1:round(80/h)));
CR_mean_inh_long=mean(CR(S_L==2,1:round(80/h)));

save('CondInhSim','CR_mean_exc','CR_mean_inh_short','CR_mean_inh_long','h')

figure
plot(h:h:length(CR_mean_inh_short)*h,CR_mean_inh_short,h:h:length(CR_mean_exc)*h,CR_mean_exc,'LineWidth',6)
PlotProperties
legend({'excitor & short inhibitor','excitor'},'Box','off')
title('model')
xlabel('time (sec)')
ylabel('response strength')

figure
plot(h:h:length(CR_mean_inh_long)*h,CR_mean_inh_long,h:h:length(CR_mean_exc)*h,CR_mean_exc,'LineWidth',6)
PlotProperties
legend({'excitor & long inhibitor','excitor'},'Box','off')
title('model')
xlabel('time (sec)')
ylabel('response strength')




