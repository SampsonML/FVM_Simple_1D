function [phi_fvm,ErrorTracker] = FVM_Func(L,D,u,c0,N,theta,sigma,lambda,BCO,BCL,dt,tend,step,plt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Matt Sampson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
%% Making the mesh
% Spatial mesh parameters for solution evaluation
x=linspace(0,L,N+1)';

% Set FVM mesh geometric properties
Dx = zeros(N,1); % distance between nodes.
dx = zeros(N,1); % Control Volume size.
Dx(1)   =  x(2)-x(1);
dx(1) = (x(2)-x(1))/2; 
for i = 2:N
    Dx(i) = x(i+1) - x(i);
    dx(i) = (x(i+1) - x(i-1))/2;
end
dx(N+1) = (x(N+1)-x(N))/2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Making Boundary Conditions
% Generalising boundary conditions
% Dirichelet Boundary x = 0
A0 = BCO(1);
B0 = BCO(2);
C0 = BCO(1)*BCO(3);

% Generalised Boundary x = L
AL = BCL(1);
BL = BCL(2);
CL = BCL(1)*BCL(3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Creating Time Steps and Exact Solution
% Time steps for plotting
% tm=[0; 0.5; 1.0; 1.5;2.0;2.5;3.0];
tm=linspace(0,tend,(tend/step + 1))';
mt=size(tm,1);

% Number of time-steps for FVM
tfinal=tm(end);
no_time_steps = tfinal/dt;


phi_exact=zeros(N+1,mt);
% Define BC for Exact
phi_exact(1,:) = 1;
phi_fvm = phi_exact;
% Forming phi for FVM
phi = phi_fvm(:,1);
ip =2;

% Loop over time to generate the derived analytic solution
for k=2:mt
        

    phi_exact(:,k)= ANALYTIC_ADVECTDIFF(x,tm(k),u,D,lambda,c0);
    
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Finite Volume Section
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Form our FVM Matrices

A = sparse(N+1,N+1);  % Implicit terms N+1
B = zeros(N+1,1);  % Explicit terms N


%% Loop Type 2


% Process time stepping for FVM and generate approximate solution
for j = 1:no_time_steps
    
    % Generalized boundary conditions x = 0
   termP = (-u*sigma/2 + D*A0/B0 + + D/Dx(2) + lambda*dx(1)); % Phi P
   termE = (-u*sigma/2 + D/Dx(2)); % Phi E
   
   fracNP = dt*theta/dx(1);
   fracN  = dt*(1-theta)/dx(1);

   A(1,1) =  1 + fracNP*termP;
   A(1,2) =  -fracNP*termE;

   B(1) = phi(1)*(1 - fracN*termP) + ...
       phi(2)*(fracN*termE) + (dt*D/dx(1))*(C0/B0);
    
    % internal control volumes dx = deltax Dx = Deltax
    for p = 2:N
        
      % Evaluate your finite volume matrix here. This matrix is called A.
       termP =  (u*(1-sigma)+ D/Dx(p-1) + D/Dx(p) + lambda*dx(p)); % Phi P
       termW =  (u*(1- sigma/2) + D/Dx(p-1)); % Phi W
       termE =   (D/Dx(p) - u*sigma/2); % Phi E
       
       fracNP = dt*theta/dx(p);
       fracN  = dt*(1-theta)/dx(p);
       
        % LHS Matrix A
        A(p,p-1) = -fracNP*termW;  % West
        A(p,p) = 1 + fracNP*termP ; % Central
        A(p,p+1) = -fracNP*termE;  % East
        
        
        % RHS Vector B
        B(p) = phi(p)*(1 - fracN*termP)+ phi(p-1)*(fracN*termW) + ...
            phi(p+1)*(fracN*termE);
        
    end
    
    % Generalised conditions at right boundary
    termP = u*(1 - sigma/2)+ D/Dx(N)+ D*AL/BL + dx(N+1)*lambda ;
    termW = u*(1 - sigma/2) + D/Dx(N);
    
    fracNP = dt*theta/dx(N+1);
    fracN =dt*(1-theta)/dx(N+1);

    A(N+1,N+1) = 1 + fracNP*termP ;
    A(N+1, N) =  -fracNP*termW;

    B(N+1) =  phi(N+1)*(1 - fracN*termP) +...
        phi(N)*(fracN*termW) +...
        (D*CL/BL)*dt/dx(N+1);
    
    
    phi=A\B;  % Solve linear system to obtain fvm solution
    
    
    % Check to see if we need to plot solution
    if (abs(j*dt - tm(ip)) < 1e-5)
        err=norm(phi-phi_exact(:,ip),inf);
        fprintf('Solution computed at t= %g error=%g \n',j*dt, err);
        phi_fvm(:,ip)=phi;
        ip=ip+1;
    end
    
end
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot solution distributions

if plt == 1

 figure;
 col = lines(mt); % Colors
 p1 = zeros(mt,1);
 for i = 1:mt
    p1((2*i-1)) = plot(x,phi_exact(:,i),'o' ,'Color',col(i,:),'LineWidth',2,'MarkerSize',4);
    hold on
    p1(2*i) = plot(x,phi_fvm(:,i),'-' ,'Color',col(i,:),'LineWidth',2,'MarkerSize',4);
    hold on
 end
% Legend displaying times
 labels = cell(mt,1);
 for i = 1:mt
    labels{2*i -1} = ['Exact Sol at t = ',num2str(tm(i))];
    labels{2*i} = ['FVM Approx at t = ',num2str(tm(i))];
 end
 legend(p1,labels,'Location','southwest','Interpreter','latex','FontSize',22);
 grid on
 grid minor
 hold off;

 xlabel('Distance from Source (Scaled Units)','FontSize',27,'FontAngle','italic','Interpreter','latex'); xlim([0 L]);
 ylabel('Concerntration of Contaminent c(x,t)','FontSize',27,'FontAngle','italic','Interpreter','latex');  ylim([0,1]);

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

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Error Plots

ErrorTracker = zeros(mt,1);

for m = 1:mt
  ErrorTracker(m) = norm(phi_fvm(:,m)-phi_exact(:,m),inf);
end

if plt ==1
 figure;
 plot(tm,ErrorTracker,'LineWidth',2)
 hold on
 plot(tm,ErrorTracker,'o','color','r','Markersize',6)
 grid on
 grid minor
 hold off;

 xlabel('Time (Days)','FontSize',27,'FontAngle','italic','Interpreter','latex'); xlim([0 tend]);
 ylabel('Error of FVM Solution $\infty$ norm','FontSize',27,'FontAngle','italic','Interpreter','latex');

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


 title({[lab1,' - ',lab2, ' Error']
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
 hold off
end
end
