% Filename: PDE_solver.m
% Author: Maud El-Hachem
% Queensland University of Technology, Brisbane, Australia, July 2021
% Reference:  M. El-Hachem, S.W. McCue, M.J. Simpson (2021) 
% Non-vanishing sharp-fronted travelling wave solutions of the
% Fisher-Kolmogorov model.
% The script contains:
%   - two calls to the function FisherStefanNonUniformGrid to generate 
%     Figures 2(a) and 2(h).
%	- the function FisherStefanNonUniformGrid
%   - the function tridia
% This function generates the time-dependent solutions
% to Equations (2.1)-(2.3) using parameters kappa, uf, initial length 
% and beta (a parameter used in the initial condition).
% A geometric, non uniform, grid is used to discretise
% the x-domain. The function plots  time-dependent solutions
% at t=5,10,15 and 20.
% The function returns the speed of the moving boundary, this an estimation
% of the wave speed, once  time-dependent solutions has reached constant
% shape and speed.

% Generating Figure 2(a)
% Calling the function FisherStefanNonUniformGrid(kappa, uf, s0, beta) 
% with the following input: kappa = 30.910, uf = 0.75, s0 = 1, beta = 0
% c1 is the estimated speed of the moving boundary
c1 = FisherStefanNonUniformGrid(30.910, 0.75, 1, 0);

% Generating Figure 2(h)
% Calling the function FisherStefanNonUniformGrid(kappa, uf, s0, beta) 
% with the following input: kappa = -2.579, uf = 0.25, s0 = 200, beta = 195
% c2 is the estimated speed of the moving boundary
c2 = FisherStefanNonUniformGrid(-2.579, 0.25, 200, 195);


% Function FisherStefanNonUniformGrid
% This function solves Equations (2.1)-(2.3) using Newton-Raphson method
% with a forward finite difference in time and a non uniform grid for x.
% The function displays the initial condition in red and  time-dependent
% solutions at t=5,10,15 and 20, and shows parameters kappa and 
% uf used in boundary condition at s(t). 
% The function estimated and displays the speed of the moving boundary.
% OUPUT ARGUMENTS: c, estimated wave speed 
% ** c = (s(T)-s(T-1))/dt, where T is the end time.
% INPUT ARGUMENTS:
% ** kappa, parameter in the Stefan condition, Equation (2.2)
% ** uf,    parameter in the boundary condition at s(t), Equation (2.2)
% ** s0,    initial length of the domain, Equation (2.3)
% ** beta,  parameter in the initial condition, Equation (2.3)
function c = FisherStefanNonUniformGrid(kappa, uf, s0, beta)

    % color used to display density profiles at t=5,10,15,20   
    ColorBlue = [0, 0.4470, 0.7410];
    N=5000; % N is the number of intervals in the non uniform grid
    dxi0 = 1e-6; % smallest grid space on RHS
    % array of position of grid nodes, 
    % that will be on a scale between 0 and 1
    xi=zeros(1,N+1);
    % array of spacings
    dxi=zeros(1,N);   
    % the last node is at xi=1
    xi(1,N+1)=1.0;
    % creating a non uniform grid, between 0 and 1,
    % a geometric grid, where the ratio of the distance beween
    % two consective spacings is constant
    % the smallest grid space must be equal to dxi0 and must be 
    % on the right hand side of the grid 
    ff = @(x) (1-x^N)/(1-x)-1/dxi0;
    ee = fsolve(ff,1.01);
    term=0;
    for i=N:-1:1
        term=term+1;
        xi(1,i) = xi(1,i+1)-dxi0*ee^(term-1);
        dxi(1,i) = xi(1,i+1)-xi(1,i);
    end
    % the first node is at xi=0
    xi(1,1)=0;  

    % dt, length of interval between two time steps
    dt=1e-3;
    % T, maximum time
    T=20;
    % current time
    t=0.0;
    % tolerance factor used in Newton-Raphson algorithm
    tol=1e-10;
    % maximum number of time steps
    maxsteps=round(T/dt);
    
    % current length of the domain, rescaled position of  moving boundary,
    % initialised at the initial length s0
    l=s0;
    % previous length, initialised at the initial length s0
    pl=l;

    % an array to record the current length at each time step
    Lrecord=zeros(maxsteps+1,1);
    Lrecord(1,1)=s0;
    % an array to record the current speed of the moving boundary
    % at each time step
    wsrecord=zeros(maxsteps+1,1);
    wsrecord(1,1)=0;
    

    % creation of matrix of coefficients used in Newton-Raphson method
    a=zeros(1,N+1);
    b=zeros(1,N+1);
    c=zeros(1,N+1);
    d=zeros(1,N+1);
    delu=ones(1,N+1);
    
    % array of density u(x,t) for each position x at current time 
    u=ones(1,N+1);
    % array of density u(x,t) for each position x at previous time
    pu=ones(1,N+1);
    
    % intialisation of density profile at t=0
    % using Equation (2.3) and parameters uf, beta and s0
    slope = -(1-uf)/(s0-beta);
    for i=1:N
        if (xi(1,i)*l)<(beta)         
            u(1,i)=1;
            pu(1,i) = 1;
        else
            u(1,i)=slope*xi(1,i)*l+(uf-slope*l);
            pu(1,i)=slope*xi(1,i)*l+(uf-slope*l); 
        end
    end
    % setting density profile to uf at boundary condition s(t)
    u(1,N+1)= uf;
    pu(1,N+1)= uf;

    figure
    hold on
    % displaying the initial condition in red
    plot(xi*l,u,'r','LineWidth',2);
    
            
    % Newton-Raphson algorithm -- main loop
    % for all time steps
    for i=1:maxsteps
        t=t+dt;
        kk=0;
        delu=ones(1,N+1);

        while norm(delu,Inf) > tol
            kk=kk+1;

            % boundary conditions at x=0, Equation (2.2)
            a(1,1)=0.0;
            b(1,1)=-1.0;
            c(1,1)=1.0;
            d(1,1)=-1*(u(1,2)-u(1,1));

            % boundary conditions at x=s(t), Equation (2.2)
            a(1,N+1)=0.0;
            b(1,N+1)=1.0;
            c(1,N+1)=0.0;
            d(1,N+1)=-(u(1,N+1)-uf);
            
            for j=2:N
                hp=dxi(j);
                hm=dxi(j-1);
                alpha = 1/(hm*(hm+hp));
                beta = -1/(hm*hp);
                delta = 1/(hp*(hm+hp));

                a(1,j)=2*alpha/l^2-xi(1,j)*(l-pl)*alpha*hp/(l*dt);
                b(1,j)=-1.0/dt+1.0-2.0*u(1,j)+2*beta/l^2-xi(1,j)*(l-pl)*beta*(hm-hp)/(l*dt);
                c(1,j)=2*delta/l^2+xi(1,j)*(l-pl)*delta*hm/(l*dt);
                d(1,j)=(u(1,j)-pu(1,j))/dt-2*(alpha*u(1,j-1)+beta*u(1,j)+delta*u(1,j+1))/l^2 ...
                -xi(1,j)*(l-pl)*(delta*hm*u(1,j+1)+beta*(hm-hp)*u(1,j)-alpha*hp*u(1,j-1))/(l*dt)...
                -u(1,j)*(1-u(1,j)); 
            end

            delu = thomas(N+1,a,b,c,d);
            % correcting the current density profile
            u(1,:)=u(1,:)+delu(1,:); 
            
            % Stefan condition, at the moving boundary, Equation (2.2)
            hp=dxi(N);
            hm=dxi(N-1);
            alpha = 1/(hm*(hm+hp));
            beta = -1/(hm*hp);
            delta = 1/(hp*(hm+hp));
            l=pl-(dt*kappa/(pl))*(uf*hm*delta...
                +u(1,N)*(hm-hp)*beta-u(1,N-1)*hp*alpha);
        end

        % displaying the solution at requested times
        if mod(i,5/dt)==0
            fprintf('Time %d\n',t);
            fprintf('Iteration %d\n',kk);
            plot(xi*l,u,'LineWidth',2,'Color',ColorBlue);
        end
        
        % updating the array of lengths for each time
        Lrecord(i+1,1)=l; 
        % updating the array of moving boudary speeds
        wsrecord(i+1,1)=(l-pl)/dt;
        % updating the current solution
        pu(1,:)=u(1,:);
        % updating the current length of the domain
        pl=l;

    end
    
    % moving boundary speed at end time
    c = wsrecord(end);
    
    % displaying important informations on the figure
    % c, kappa and uf
    textc = strcat(strcat('$c = ',num2str(c,3)), '$');
    text(50,0.1,textc,'interpreter','latex','fontsize',24)
    textc = strcat(strcat('$\kappa = ',num2str(kappa,3)), '$');
    text(50,0.3,textc,'interpreter','latex','fontsize',24)
    textc = strcat(strcat('$u_\textrm{f} = ',num2str(uf,3)), '$');
    text(50,0.5,textc,'interpreter','latex','fontsize',24)
    % displaying x and y labels
    xlabel('$x$','interpreter','latex','fontsize',24)
    ylabel('$u(x,t)$','interpreter','latex','fontsize',24)
    % enhancing the figure
    hold off
    box on
    set(groot, 'defaultAxesTickLabelInterpreter','latex');
    set(groot, 'defaultLegendInterpreter','latex');
    set(gca,'fontsize', 24);
    ylim([0 1])
    xlim([0 200])
end

%% Function tridia
% This function implements Thomas algorithm that solves a system
% of equations Ax = d, where A is a tridiagonal matrix. The parameters 
% a,b and c are three diagonals of matrix A. N is size of 
% vector solution x.
function x = thomas(N,a,b,c,d)
x=zeros(1,N);
    bb=b;
    dd=d;
    for i=2:N
        ff=a(i)/bb(i-1);
        bb(i)=bb(i)-c(i-1)*ff;
        dd(i)=dd(i)-dd(i-1)*ff;
    end
    
    for i=1:N-1
    x(N)=dd(N)/bb(N);    
    j=N-i;
    x(j)=(dd(j)-c(j)*x(j+1))/bb(j);
    end
end



