%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MXB326 Assignment 1
% Run Script Matt Sampson
% Produces all plots used in report (and a few more the were decided not to
% use.

%% Section 3 Computational Model --- 2a Homogeneous Direchlet x = L)
L = 1.5;          % Right boundary scaled length
D = 0.002;        % Diffusion constant
u = 0.4;          % Flow speed in stream
c0 = 1;           % Concerntration released at factory site
N = 300;          % Points in Mesh
thetaBE = 1;      % 1 = Back Euler 
thetaCN = 1/2;    % 1/2 = Crank Nicolson
sigmaAV = 1;      % 1 = Averaging
sigmaUW = 0;      % 0 = Upwinding
lambda = 0.05;    % Source term
BCO = [1e7,1,1];  % Boundary Conditions x = 0
BCL = [1e7,1,0];  % Boundary Conditions x = L Homogenous Dirichlet
dt = 2e-4;        % Time step length
tend = 4;         % Final time to compute solutions at
step = 0.5;       % Time step for calculations
plt = 1;          % Plots turned to true


% Trial 1 BE-UW
[phiBU,ErrorBU] = FVM_Func(L,D,u,c0,N,thetaBE,sigmaUW,lambda,BCO,BCL,dt,tend,step,plt);

% Trial 2 BE-AV
[phiBA, ErrorBA] = FVM_Func(L,D,u,c0,N,thetaBE,sigmaAV,lambda,BCO,BCL,dt,tend,step,plt);

% Trial 3 CN-UW
[phiCU,ErrorCU] = FVM_Func(L,D,u,c0,N,thetaCN,sigmaUW,lambda,BCO,BCL,dt,tend,step,plt);

% Trial 4 CN-AV
[phiCA,ErrorCA] = FVM_Func(L,D,u,c0,N,thetaCN,sigmaAV,lambda,BCO,BCL,dt,tend,step,plt);

%% Plotting the error in one plot
plt = 0; % Turn plots off, make FVM data with smaller timesteps
step = 0.05;
[phiBU,ErrorBU] = FVM_Func(L,D,u,c0,N,thetaBE,sigmaUW,lambda,BCO,BCL,dt,tend,step,plt);
[phiBA, ErrorBA] = FVM_Func(L,D,u,c0,N,thetaBE,sigmaAV,lambda,BCO,BCL,dt,tend,step,plt);
[phiCU,ErrorCU] = FVM_Func(L,D,u,c0,N,thetaCN,sigmaUW,lambda,BCO,BCL,dt,tend,step,plt);
[phiCA,ErrorCA] = FVM_Func(L,D,u,c0,N,thetaCN,sigmaAV,lambda,BCO,BCL,dt,tend,step,plt);


time = linspace(0,tend,(tend/step+1));

figure;
plot(time,ErrorBU,'LineWidth',2)
hold on
plot(time,ErrorBA,'-','color','r','LineWidth',2)
plot(time,ErrorCU,'o','color','b','LineWidth',2)
plot(time,ErrorCA,'o','color','g','LineWidth',2)
grid on
grid minor
hold off;

xlabel('Time (Days)','FontSize',27,'FontAngle','italic','Interpreter','latex');
ylabel('Error of FVM Solution $\infty$ norm','FontSize',27,'FontAngle','italic','Interpreter','latex');


title({['Error of FVM Solutions Dirichlet Bounds']
    ['D = ' num2str(D) '  L = ' num2str(L)...
    ' $\lambda$ = ' num2str(lambda)]},...
        'FontSize',29,'FontAngle','italic','Interpreter','latex');

legend('BE-UW', 'BE-AV','CN-UW','CN-AV',...
    'Location','northeastoutside','Interpreter','latex','FontSize',22);
    
% Extra Plot Aesthetics
ax = gca; % current axes
ax.FontSize = 26;
ax.TickDir = 'out';
ax.TickLength = [0.01 0.01];
ax.TickLabelInterpreter = 'latex';
ax.XMinorTick = 'on';
hold off


%% Section 3 Computational Model --- 2b Homogeneous Neumann)
step = 0.5;
BCO = [1e7,1,1];  % Boundary Conditions x = 0
BCL = [0,1,0];  % Boundary Conditions x = L Homogeneous Neumann
plt = 1;

% Trial 5 BE-UW
[phiBU,ErrorBU] = FVM_Func(L,D,u,c0,N,thetaBE,sigmaUW,lambda,BCO,BCL,dt,tend,step,plt);

% Trial 6 BE-AV
[phiBA, ErrorBA] = FVM_Func(L,D,u,c0,N,thetaBE,sigmaAV,lambda,BCO,BCL,dt,tend,step,plt);

% Trial 7 CN-UW
[phiCU,ErrorCU] = FVM_Func(L,D,u,c0,N,thetaCN,sigmaUW,lambda,BCO,BCL,dt,tend,step,plt);

% Trial 8 CN-AV
[phiCA,ErrorCA] = FVM_Func(L,D,u,c0,N,thetaCN,sigmaAV,lambda,BCO,BCL,dt,tend,step,plt);


%% Plotting the error on one plot Neumann Bounds
plt = 0; % Turn plots off, make FVM data with smaller timesteps
step = 0.05;
[phiBU,ErrorBU] = FVM_Func(L,D,u,c0,N,thetaBE,sigmaUW,lambda,BCO,BCL,dt,tend,step,plt);
[phiBA, ErrorBA] = FVM_Func(L,D,u,c0,N,thetaBE,sigmaAV,lambda,BCO,BCL,dt,tend,step,plt);
[phiCU,ErrorCU] = FVM_Func(L,D,u,c0,N,thetaCN,sigmaUW,lambda,BCO,BCL,dt,tend,step,plt);
[phiCA,ErrorCA] = FVM_Func(L,D,u,c0,N,thetaCN,sigmaAV,lambda,BCO,BCL,dt,tend,step,plt);

time = linspace(0,tend,(tend/step+1));

figure;
plot(time,ErrorBU,'LineWidth',2)
hold on
plot(time,ErrorBA,'-','color','r','LineWidth',2)
plot(time,ErrorCU,'o','color','b','LineWidth',2)
plot(time,ErrorCA,'o','color','g','LineWidth',2)
grid on
grid minor
hold off;

xlabel('Time (Days)','FontSize',27,'FontAngle','italic','Interpreter','latex');
ylabel('Error of FVM Solution $\infty$ norm','FontSize',27,'FontAngle','italic','Interpreter','latex');


title({['Error of FVM Solutions Neumann Bounds']
    ['D = ' num2str(D) '  L = ' num2str(L)...
    ' $\lambda$ = ' num2str(lambda)]},...
        'FontSize',29,'FontAngle','italic','Interpreter','latex');

legend('BE-UW', 'BE-AV','CN-UW','CN-AV',...
    'Location','northeastoutside','Interpreter','latex','FontSize',22);
    
% Extra Plot Aesthetics
ax = gca; % current axes
ax.FontSize = 26;
ax.TickDir = 'out';
ax.TickLength = [0.01 0.01];
ax.TickLabelInterpreter = 'latex';
ax.XMinorTick = 'on';
hold off


%% Section 4 Analysis

% Question (i) Lethal Dose
% Numerical Integration
tend = 4;
step = 0.025;
[Phi_FVM] = FVM_Func(L,D,u,c0,N,thetaCN,sigmaAV,lambda,BCO,BCL,dt,tend,step,plt);

Integral = zeros(tend/step,1);

for j = 2:(tend/step)
      Integral(j) = Integral(j-1) + step*Phi_FVM(200,j); % Add c(1,t) part
end

for i = 1:length(Integral(:,1))
    %Time = num2str(i*step);
    disp(['Concerntration at t= ',num2str(i*step),' is:  ']);disp(Integral(i));
end


%% Plotting Integral

time1 = linspace(0,tend,(tend/step+1));
time2 = linspace(0,tend,(tend/step)); % Calculate one less step in integral
figure;
subplot(1,2,1)
plot(time1,Phi_FVM(200,:),'LineWidth',2)
hold on
grid on
grid minor
xlabel('Time (Days)','FontSize',27,...
    'FontAngle','italic','Interpreter','latex');
ylabel('Concerntration Level ',...
    'FontSize',27,'FontAngle','italic','Interpreter','latex');
title({['Concerntration Level $x = 1$']
    ['D = ' num2str(D) '  L = ' num2str(L)...
    ' $\lambda$ = ' num2str(lambda)]},...
        'FontSize',29,'FontAngle','italic','Interpreter','latex');

subplot(1,2,2)
plot(time2,Integral,'LineWidth',2)
hold on
yline(0.8,'LineWidth',2,'color','r')
grid on
grid minor
xlabel('Time (Days)','FontSize',27,...
    'FontAngle','italic','Interpreter','latex');
ylabel('Concerntration Level',...
    'FontSize',27,'FontAngle','italic','Interpreter','latex');

title({['Cumulative Concerntration Level $x = 1$']
    ['D = ' num2str(D) '  L = ' num2str(L)...
    ' $\lambda$ = ' num2str(lambda)]},...
        'FontSize',29,'FontAngle','italic','Interpreter','latex');
    
% % Extra Plot Aesthetics
% ax = gca; % current axes
% ax.FontSize = 26;
% ax.TickDir = 'out';
% ax.TickLength = [0.01 0.01];
% ax.TickLabelInterpreter = 'latex';
% ax.XMinorTick = 'on';
% hold off

%% Question (ii)

step = 0.5;
tend = 3;
dt = 2e-4;
N = 300;
BCO = [1e5,1,1];  % Boundary Conditions x = 0
BCL = [0,1,0];  % Boundary Conditions x = L Homogeneous Neumann

[phi_test] = FVM_Test(L,D,u,c0,N,thetaCN,sigmaAV,lambda,BCO,BCL,dt,tend,step,plt);

[phi_com] = FVM_Func(L,D,u,c0,N,thetaCN,sigmaAV,lambda,BCO,BCL,dt,tend,step,plt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot 2 solution distributions

tm=linspace(0,tend,(tend/step + 1))';
mt=size(tm,1);
x=linspace(0,L,N+1)';
figure;
col = lines(mt); % Colors

p1 = zeros(mt,1);
for i = 1:mt
    
    p1((2*i-1)) = plot(x,phi_test(:,i),':' ,'Color',col(i,:),'LineWidth',2,'MarkerSize',4);
    hold on
    p1(2*i) = plot(x,phi_com(:,i),'-' ,'Color',col(i,:),'LineWidth',2,'MarkerSize',4);
    hold on
   
end
% Legend displaying times
labels = cell(mt,1);
for i = 1:mt
    labels{2*i -1} = ['c(0,t) = $\exp\{-0.9(t - 0.25) \} \ t > 0.25$: t = ',num2str(tm(i))];
    labels{2*i} = ['c(0,t) = 1: t = ',num2str(tm(i))];
end
legend(p1,labels,'Location','northeastoutside','Interpreter','latex','FontSize',22);
grid on
grid minor
hold off;

xlabel('Distance from Source (Scaled Units)','FontSize',27,'FontAngle','italic','Interpreter','latex'); xlim([0 L]);
ylabel('Concerntration of Contaminent c(x,t)','FontSize',27,'FontAngle','italic','Interpreter','latex');  ylim([0,1]);


theta = 1/2;
sigma = 1;
% Title auto-generate title
if theta == 1
    lab1 = 'Back Euler';
elseif theta == 1/2
    lab1 = 'Crank Nicolson';
end
if sigma == 1
    lab2 = 'Averaging';
elseif sigma == 0
    lab2 = 'Upwinding';
end


title({[lab1,' - ',lab2, ' Advection Model']
    ['D = ' num2str(D) '  L = ' num2str(L)...
    ' $\lambda$ = ' num2str(lambda)]},...
        'FontSize',29,'FontAngle','italic','Interpreter','latex');
    
% Extra Plot Aesthetics
ax = gca; % current axes
ax.FontSize = 26;
ax.TickDir = 'out';
ax.TickLength = [0.01 0.01];
ax.TickLabelInterpreter = 'latex';
ax.XMinorTick = 'on';

%% Integral (ii) 
step = 0.05;
tend = 4;

[phi_test] = FVM_Test(L,D,u,c0,N,thetaCN,sigmaAV,lambda,BCO,BCL,dt,tend,step,plt);
Integral = zeros(tend/step,1);

for j = 2:(tend/step)
      Integral(j) = Integral(j-1) + step*phi_test(200,j); % Add c(1,t) part
end

for i = 1:length(Integral(:,1))
    %Time = num2str(i*step);
    disp(['Concerntration at t= ',num2str(i*step),' is:  ']);disp(Integral(i));
end

%%%%%%
%%
time2 = linspace(0,tend,(tend/step)); % Calculate one less step in integral
figure
plot(time2,Integral,'LineWidth',2)
hold on
yline(0.8,'LineWidth',2,'color','r')
grid on
grid minor
xlabel('Time (Days)','FontSize',27,...
    'FontAngle','italic','Interpreter','latex');
ylabel('Concerntration Level',...
    'FontSize',27,'FontAngle','italic','Interpreter','latex');

title({['Cumulative Concerntration Level $x = 1$']
    ['D = ' num2str(D) '  L = ' num2str(L)...
    ' $\lambda$ = ' num2str(lambda)]},...
        'FontSize',29,'FontAngle','italic','Interpreter','latex');

%% Question (iii)

step = 0.05;
tend = 5;
lambda = 0.3;

BCO = [1e7,1,1];  % Boundary Conditions x = 0
BCL = [0,1,0];  % Boundary Conditions x = L Homogeneous Neumann

% Numerical Integration
[Phi_FVM] = FVM_Func(L,D,u,c0,N,thetaCN,sigmaAV,lambda,BCO,BCL,dt,tend,step,plt);

Integral = zeros(tend/step,1);

for j = 2:(tend/step)
      Integral(j) = Integral(j-1) + step*Phi_FVM(200,j); % Add c(1,t) part
end

for i = 1:length(Integral(:,1))
    %Time = num2str(i*step);
    disp(['Concerntration at t= ',num2str(i*step),' is:  ']);disp(Integral(i));
end


%% Plotting Integral
time2 = linspace(0,tend,(tend/step)); % Calculate one less step in integral
figure;
plot(time2,Integral,'LineWidth',2)
hold on
yline(0.8,'LineWidth',2,'color','r')
grid on
grid minor
xlabel('Time (Days)','FontSize',27,...
    'FontAngle','italic','Interpreter','latex');
ylabel('Concerntration Level',...
    'FontSize',27,'FontAngle','italic','Interpreter','latex');

title({['Cumulative Concerntration Level $x = 1$']
    ['D = ' num2str(D) '  L = ' num2str(L)...
    ' $\lambda$ = ' num2str(lambda)]},...
        'FontSize',29,'FontAngle','italic','Interpreter','latex');
    
%% Q3 part two, new lambda time dependant boundary

lambda = 0.3;
step = 0.5;
tend = 4;
dt = 2e-4;
N = 300;
BCO = [1e5,1,1];  % Boundary Conditions x = 0
BCL = [0,1,0];  % Boundary Conditions x = L Homogeneous Neumann

[phi_test] = FVM_Test(L,D,u,c0,N,thetaCN,sigmaAV,lambda,BCO,BCL,dt,tend,step,plt);

[phi_com] = FVM_Func(L,D,u,c0,N,thetaCN,sigmaAV,lambda,BCO,BCL,dt,tend,step,plt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot 2 solution distributions

tm=linspace(0,tend,(tend/step + 1))';
mt=size(tm,1);
x=linspace(0,L,N+1)';
figure;
col = lines(mt); % Colors

p1 = zeros(mt,1);
for i = 1:mt
    
    p1((2*i-1)) = plot(x,phi_test(:,i),':' ,'Color',col(i,:),'LineWidth',2,'MarkerSize',4);
    hold on
    p1(2*i) = plot(x,phi_com(:,i),'-' ,'Color',col(i,:),'LineWidth',2,'MarkerSize',4);
    hold on
   
end
% Legend displaying times
labels = cell(mt,1);
for i = 1:mt
    labels{2*i -1} = ['c(0,t) = $\exp\{-0.9(t - 0.25) \} \ t > 0.25$: t = ',num2str(tm(i))];
    labels{2*i} = ['c(0,t) = 1: t = ',num2str(tm(i))];
end
legend(p1,labels,'Location','northeastoutside','Interpreter','latex','FontSize',22);
grid on
grid minor
hold off;

xlabel('Distance from Source (Scaled Units)','FontSize',27,'FontAngle','italic','Interpreter','latex'); xlim([0 L]);
ylabel('Concerntration of Contaminent c(x,t)','FontSize',27,'FontAngle','italic','Interpreter','latex');  ylim([0,1]);


theta = 1/2;
sigma = 1;
% Title auto-generate title
if theta == 1
    lab1 = 'Back Euler';
elseif theta == 1/2
    lab1 = 'Crank Nicolson';
end
if sigma == 1
    lab2 = 'Averaging';
elseif sigma == 0
    lab2 = 'Upwinding';
end


title({[lab1,' - ',lab2, ' Advection Model']
    ['D = ' num2str(D) '  L = ' num2str(L)...
    ' $\lambda$ = ' num2str(lambda)]},...
        'FontSize',29,'FontAngle','italic','Interpreter','latex');
    
% Extra Plot Aesthetics
ax = gca; % current axes
ax.FontSize = 26;
ax.TickDir = 'out';
ax.TickLength = [0.01 0.01];
ax.TickLabelInterpreter = 'latex';
ax.XMinorTick = 'on';

%% Integral (iii) 
step = 0.05;
tend = 10;

[phi_test] = FVM_Test(L,D,u,c0,N,thetaCN,sigmaAV,lambda,BCO,BCL,dt,tend,step,plt);
Integral = zeros(tend/step,1);

for j = 2:(tend/step)
      Integral(j) = Integral(j-1) + step*phi_test(200,j); % Add c(1,t) part
end

for i = 1:length(Integral(:,1))
    %Time = num2str(i*step);
    disp(['Concerntration at t= ',num2str(i*step),' is:  ']);disp(Integral(i));
end

%%%%%%
%%
time2 = linspace(0,tend,(tend/step)); % Calculate one less step in integral
figure
plot(time2,Integral,'LineWidth',2)
hold on
yline(0.8,'LineWidth',2,'color','r')
grid on
grid minor
xlabel('Time (Days)','FontSize',27,...
    'FontAngle','italic','Interpreter','latex');
ylabel('Concerntration Level',...
    'FontSize',27,'FontAngle','italic','Interpreter','latex');

title({['Cumulative Concerntration Level $x = 1$']
    ['D = ' num2str(D) '  L = ' num2str(L)...
    ' $\lambda$ = ' num2str(lambda)]},...
        'FontSize',29,'FontAngle','italic','Interpreter','latex');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of the run script
