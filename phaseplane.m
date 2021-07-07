% Filename: phasePlane.m
% Author: Maud El-Hachem
% Queensland University of Technology, Brisbane, Australia, July 2021
% Reference:  M. El-Hachem, S.W. McCue, M.J. Simpson (2021) 
% Non-vanishing sharp-fronted travelling wave solutions of the
% Fisher-Kolmogorov model.
% The script contains:
%   - two calls to the function plotPhasePlane to generate 
%     Figures 3(a) and (g).
%	- the function plotPhasePlane
%   - the function HeunsSolver

% Generating Figure 3(a)
% Calling the function plotPhasePlane(c,Uf) where c = 2.5 and Uf = 0.5
% kappa1_num is the numerical kappa that we obtain from
% the Stefan condition in the plane, such as kappa = -2.5/V(U==0.5)
kappa1_num = plotPhasePlane(2.5, 0.5);

% Generating Figure 3(g)
% Calling the function plotPhasePlane(c,Uf) where c = -1 and Uf = 0.5
kappa2_num = plotPhasePlane(-1, 0.5);
% kappa2_num is the numerical kappa that we obtain from
% the Stefan condition in the plane, such as kappa = 1/V(U==0.5)

% Function plotPhasePlane
% This function solves Equations (2.6) and (2.7) in the phase plane with
% Heun's method and plots the solution on the plane V(z) versus U(z).
% It also show the solution trunctated at the intersection point where
% c = -kappa * V(U==Uf).
% The same plot shows the equilibrium points (0,0) and (1,0) 
% and the intersection point of the solution with the axis U(z)=Uf.
% OUTPUT ARGUMENTS: kappa_num, the numerical kappa that we obtain from
% the Stefan condition in the plane c = -kappa * V(U==Uf)
% INPUT ARGUMENTS: c, the wave speed 
function kappa_num = plotPhasePlane(c,Uf)

    % colours used to display the solutions
    ColorOrange = [217 118 0]/255;
    ColorBlue = [0, 0.4470, 0.7410];
    % step size dz, z domain goes from z_begin to z_end.
    dz = 0.0001;
    z_begin = 0;
    z_end = 200;

    % depending on c, 
    % - calculating the eigenvalues and the eigenvectors of the solution 
    %   around the equilibrium point (1,0)
    % - determining the initial conditions close to the equilibrium point
    %   (1,0) along the eigenvector of the solution
    Us = 1;
    A = [0 1;(2*Us-1) -c];

    [v,d]=eig(A);
    IC = [0;0];
    if (c>0)
        if (d(1,1) > 0)
            IC = v(:,1);
        end
        if (d(2,2) > 0)
            IC = v(:,2);
        end
    end
    if (c<0)
        if (d(1,1) < 0)
            IC = v(:,1);
        end
        if (d(2,2) < 0)
            IC = v(:,2);
        end
    end

    % setting the initial conditions
    IC2 = IC(1,1)*0.0001;
    IC1 = Us+IC(2,1)*0.0001;

    % solving Equations Equations (2.6) and (2.7) in the phase plane with
    % Heun's method
    [U, V] = heunSolver(c, dz, z_begin, z_end, IC1, IC2);

    % finding the intersection point with the axis U(z)=Uf
    intersects = [];
    for ii = 1:length(U)-1
        if (U(ii,1)>= Uf && U(ii+1,1)<= Uf || U(ii,1)<= Uf && U(ii+1,1)>= Uf )
            var1 = [U(ii,1) U(ii+1,1)];
            var2 = [V(ii,1) V(ii+1,1)];
            intersects = [intersects interp1(var1,var2,Uf,'linear')];
        end
    end

    % finding V where U=Uf;
    intersect = min(intersects);    
    % calculating the numerical kappa that could be used to solve the PDE
    % kappa = -c/V(U==Uf)
    kappa_num = -c/intersect;
   
    % truncating the solution to respect the boundary condition at U(z)=Uf
    jj = 1;
    for ii = length(U):-1:1
        if (U(ii,1) > Uf && U(ii,1) < 1 && V(ii,1) < 0 && V(ii,1) > intersect)
            Utrunc(jj,1) = U(ii,1);
            Vtrunc(jj,1) = V(ii,1);
            jj = jj + 1;
        end
    end
    
    figure
    % displaying the axis U(z) and V(z), and the axis U(z)=Uf
    line([-0.25 1.1],[0 0],'Color','k','LineStyle','-','LineWidth',1);
    hold on
    line([0 0],[-5.5 0.5],'Color','k','LineStyle','-','LineWidth',1);
    line([Uf Uf],[-5.5 0.5],'Color','m','LineStyle','-','LineWidth',1);
    % displaying the sumperimposed solutions 
    % (complete trajectory and truncated trajectory)
    plot(U,V,'--','LineWidth',2,'Color',ColorOrange);
    plot(Utrunc,Vtrunc,'r-','LineWidth',2,'Color',ColorBlue);
    %displaying the intersection point with U(z) = Uf
    plot(Uf,intersect,'mo','MarkerEdgeColor','m','MarkerFaceColor','m','LineWidth',4);
    % displaying the equilibrium points (0,0) and (1,0)
    plot(0,0,'o','MarkerEdgeColor','k','MarkerFaceColor','k','LineWidth',4);
    plot(1,0,'o','MarkerEdgeColor','k','MarkerFaceColor','k','LineWidth',4);
    % displaying the value of the wave speed c
    textc = strcat(strcat('$c = ',num2str(round(c,2))), '$');
    text(0.6,intersect-0.1,textc,'interpreter','latex','fontsize',24)
    % displaying the labels of the corresponding axis
    ylabel('$V(z)$','interpreter','latex');
    xlabel('$U(z)$','interpreter','latex'); 
    % enhancing the figure
    set(groot, 'defaultAxesTickLabelInterpreter','latex');
    set(groot, 'defaultLegendInterpreter','latex');
    set(gca,'fontsize', 24);
    xlim([-0.1 1])
    ylim([intersect-0.25 0.1])
    box on;
    hold off

end

% Function heunSolver
% This function solves Equations (11) and (12) by Heun's method 
% INPUT ARGUMENTS:
% ** c, the wave speed, c can any real value
% ** dz, the step size used to discretise the domain of z
% ** z_begin and z_end, the lower and upper limit of the numerical domain 
% of z such as z_begin <= z <= z_end. The initial conditions are applied at
% z = z_begin.
% ** V1, X1, the values of the initial conditions
% OUTPUT ARGUMENTS:
% ** Uout : The solution U(z)
% ** Vout : The solution V(z)
function [Uout, Vout] = heunSolver(c,dz,z_begin,z_end,U1,V1)

    % z domain
    z = z_begin:dz:z_end;
    % number of nodes in the domain
    sz = length(z);

    % initialisation 
    V = zeros(sz,1);
    U = zeros(sz,1);
    U(1) = U1;
    V(1) = V1;
    
    % for all steps in the domain
    for i = 1:sz-1
        Ubar = dz * V(i) + U(i); 
        Vbar = dz * (-c*V(i)- U(i)*(1-U(i))) + V(i);
        U(i+1) = dz/2 * (V(i)+Vbar) + U(i); 
        V(i+1) = dz/2 * ((-c*V(i) - U(i)*(1-U(i))) + (-c*Vbar - Ubar *(1-Ubar))) + V(i); 
    end

    Uout = U;
    Vout = V;
end