%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%       Developed by Daniel Burbano Lombana, 07/05/2022
%%%   for questions contact me at daniel.burbano@rutgers.edu
%
%
%    Desciption: This code generates the two-dimensional plots of an ADN
%    with two concurrent strains
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc;
clear all;
close all;

p=60;               % control parameter
npointsGrid = 20;   % number of points of the grid (2D plot)
NumIter = 5;        % number of iterations
NT = 600;           % Total number of time steps
dt = 0.5;           % time step
N = 1e4;            % number of nodes
m = 20;             % Num of links
sigma_1 = 0.5;      % latency to become infectious of strain 1
sigma_2 = 0.5;      % latency to become infectious of strain 2
mu1 = 0.2;          % recovery rate for strain 1
mu2 = 0.2;          % recovery rate for strain 2
Ni1 = 0.01;         % Percentage if initial infected nodes with strain 1;
Ni2 = 0.01;         % Percentage if initial infected nodes with strain 2;
Wtime = 7/dt;       % time period of the control strategy ( 7 days)
WDuration = 0.5;    % duration of the home-isolation period



%% create the activity distribution
y = -2.1;          % heavy tail exponent
xmin = 0.001;      % lower cutoff of the power law
xmax = 1;          % higher cutoff of the power law
xxi = activityPotential(y,xmin,xmax,N);
eta = 10;          % Activity gain
a = eta*xxi';
a(a>1/dt) = 1/dt;  % cut-off on the activities



%% Save parameter values in a structure
IniStrain = 0.5*(1/dt); % Initial time where the epidemic process starts.
al21 = 0.1;             % strain-specific re-infection probability for strain 1
al12 = 0.1;             % strain-specific re-infection probability for strain 2


Parameters.N = N;
Parameters.m = m;
Parameters.dt = dt;
Parameters.T = NT;
Parameters.a = a;
Parameters.mu1 = mu1;
Parameters.mu2 = mu2;
Parameters.sigma_1 = sigma_1;
Parameters.sigma_2 = sigma_2;
Parameters.al12 = al12;
Parameters.al21 = al21;
Parameters.Ni1 = Ni1;
Parameters.Ni2 = Ni2;
Parameters.IniStrain = IniStrain;

Parameters.p = p;
Parameters.Wtime = Wtime;
Parameters.WDuration = WDuration;

I1 = zeros(NumIter,NT);
Inf_1 = zeros(NumIter,NT);
Inf_2 = zeros(NumIter,NT);
NssI = zeros(NumIter,NT);
Rec_1 = zeros(NumIter,NT);
Rec_2 = zeros(NumIter,NT);
Exp_1 = zeros(NumIter,NT);
Exp_2 = zeros(NumIter,NT);

PeakInf_s1 = [];
PeakInf_s2 = [];
PeakInf_tot = [];
SS_s1 = [];
SS_s2 = [];

Par1 = linspace(0,1,npointsGrid);
Par2 = linspace(0,0.5,npointsGrid);

for p1=1:length(Par1)

    Parameters.al11 = Par1(p1); % strain-specific re-infection probability for strain 1
    Parameters.al22 = Par1(p1); % strain-specific re-infection probability for strain 2

    % Display progress
    if ~mod(p1,10)
        disp(p2)
    end

    for p2=1:length(Par2)

        Parameters.lambda1 = Par2(p2);
        Parameters.lambda2 = Par2(p2);


        parfor k=1:NumIter

            out = ADN_2Var_Controlled(Parameters);

            Ns = out.Ns;           % Number of suceptible nodes
            Ne = out.Ne;           % Number of exposed nodes
            Ninf = out.Ni;         % Number of infected nodes
            Nrec = out.Nrec;       % Number of recovered
            NreExp = out.NreExp;  % number of nodes that are reexposed
            NreInf = out.NreInf;  % number of nodes that are reinfected
            NreRec = out.NreRec;  % number of nodes that are re-recovered

            Sus(k,:) = Ns;
            %
            Inf_1(k,:) = Ninf(:,1);
            Inf_2(k,:) = Ninf(:,2);
            %
            NssI(k,:) = Ns;
            %
            Rec_1(k,:) = Nrec(:,1);
            Rec_2(k,:) = Nrec(:,2);
            %
            %
            Exp_1(k,:) = Ne(:,1);
            Exp_2(k,:) = Ne(:,2);
            %
            %
            ReExp_1(k,:) = NreExp(:,1);
            ReExp_2(k,:) = NreExp(:,2);
            %
            %
            ReInf_1(k,:) = NreInf(:,1);
            ReInf_2(k,:) = NreInf(:,2);

        end

        %% Variables
        I1m = mean(Inf_1);
        I2m = mean(Inf_2);

        maxV1 = max(I1m);
        maxV2 = max(I2m);

        mV1 = mean(I1m(end-50:end));
        mV2 = mean(I2m(end-50:end));

        PeakInf_s1(p1,p2) = maxV1;
        PeakInf_s1(p1,p2) = maxV2;
        PeakInf_tot(p1,p2) = maxV1+maxV2;

        SS_s1(p1,p2) = mV1;
        SS_s2(p1,p2) = mV2;
        SS_tot(p1,p2) = mV1+mV2;

    end
end


%% Calculate stability threshold
[X,Y] = meshgrid(Par2,Par1);
%epidemic
[Lambda,Lambda2] = stability_switchedSys(Parameters,Par1,Par2);
[M,c] = contour(X,Y,Lambda,[0 0],'--w','LineWidth',3);
c.LineWidth = 2;
contourTable = getContourLineCoordinates(M);
Xc = contourTable.X;
Yc = contourTable.Y;
% endemic

[M2,c2] = contour(X,Y,Lambda2,[0 0],'--w','LineWidth',3);
c2.LineWidth = 2;
contourTable2 = getContourLineCoordinates(M2);
X2c = contourTable2.X;
Y2c = contourTable2.Y;

%% plots 
hh = 9000;
[X,Y] = meshgrid(Par2,Par1);
figure('Color', [1 1 1])
axes('Box','on', 'FontSize',25, 'FontWeight', 'normal','FontName','Arial'),hold on,
s=surf(X,Y,PeakInf_tot);
s.EdgeColor = 'none';
hold on
plot3(Xc,Yc,hh*(ones(1,length(Xc))),'--w','LineWidth',3)
title('Peak of Infected')
axis([0,0.5,0,1])
set(gca,'TickLabelInterpreter','latex')
colorbar('TickLabelInterpreter', 'latex'); 


[X,Y] = meshgrid(Par2,Par1);
figure('Color', [1 1 1])
axes('Box','on', 'FontSize',25, 'FontWeight', 'normal','FontName','Arial'),hold on,
s=surf(X,Y,SS_tot);
s.EdgeColor = 'none';
hold on
plot3(Xc,Yc,hh*(ones(1,length(Xc))),'--w','LineWidth',3)
plot3(X2c,Y2c,hh*(ones(1,length(X2c))),'--w','LineWidth',3)

title('Steady state of Infected')
axis([0,0.5,0,1])
set(gca,'TickLabelInterpreter','latex')
colorbar('TickLabelInterpreter', 'latex'); 







