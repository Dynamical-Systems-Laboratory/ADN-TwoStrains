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

function out = ADN_2Var_Controlled(Parameters)


N = Parameters.N;
dt = Parameters.dt;
T = Parameters.T;
a = Parameters.a;
ni1 = Parameters.Ni1;
ni2 = Parameters.Ni2;
IniStrain = Parameters.IniStrain;
coin = rand(N,T);


%% Initial state of infected
Nini1 = round((ni1)/100*N);
X = zeros(N,1);
Xid1 = round(1 + (N-1)*rand(Nini1,1));
X(Xid1) = 2;

Nini2 = round((ni2)/100*N);
Xid2 = round(1 + (N-1)*rand(Nini2,1));

p = Parameters.p;
Wtime = Parameters.Wtime;
WDuration = Parameters.WDuration;

IDcontrol = [];

for n=1:T    
    
    if n==IniStrain
        X(Xid2) = 2.5;
    end
    %% counting states 
    Xsu = X(X==0);                % Suceptible nodes
    Ns(n,1) = length(Xsu);        % number of suceptible nodes 
    
    Xexp1 = X(X==1);              % Exposed to 1
    Ne(n,1) = length(Xexp1);      % number of Exposed to 1
    Xexp2 = X(X==1.5);            % Exposed to 2
    Ne(n,2) = length(Xexp2);      % number of Exposed to 2
    
    Xinf1 = X(X==2);                 % Infected of 1
    Ni(n,1) = length(Xinf1);         % number of Infected of I1
    Xinf2 = X(X==2.5);               % Infected of 2
    Ni(n,2) = length(Xinf2);         % number of Infected of I2
    
    Xrec1 = X(X==3);                 % Recovered from 1
    Nrec(n,1) = length(Xrec1);       % number of Recovered from 1
    Xrec2 = X(X==3.5);               % Recovered from 2
    Nrec(n,2) = length(Xrec2);       % number of Recovered from 2
    
    Xexp12 = X(X==4);                % Exposed to strain 2 given they had strain 1 already
    NreExp(n,1) = length(Xexp12);    
    Xexp21 = X(X==4.5);              % Exposed to strain 1 given they had strain 2 already 
    NreExp(n,2) = length(Xexp21);    

    Xinf12 = X(X==5);                % Exposed to strain 2 given they had strain 1 already
    NreInf(n,1) = length(Xinf12);    
    Xinf21 = X(X==5.5);              % Exposed to strain 1 given they had strain 2 already 
    NreInf(n,2) = length(Xinf21);    

    XRecovered = X(X==6);            % Exposed to strain 2 given they had strain 1 already
    NreRec(n,1) = length(XRecovered);    
     
    Xn = X; % auxiliary variable
    XX(:,n) = X; 

    Acti = coin(:,n) <= a.*dt;


    %% Control action

    if mod(n,Wtime) == 0 || n==1
        
        nNc = round(p*N/100); % number of nodes to be controlled
        nID = randperm(N);    % Nodes to be remved
        nID(nNc+1:end) = [];
        IDcontrol = nID;
    end
    Ss(n) = (Wtime-1)*WDuration - mod(n,Wtime);
    if Ss(n) >=0
        swichingSignal(n) = 1;
        Acti(IDcontrol) = 0;
    else
        swichingSignal(n) = 0;
    end

    Xn = Node_Dynamics_Control(Acti,IDcontrol,X,Xn,Parameters);
    X = Xn;
    
end

out.Ns = Ns;          % Number of suceptible nodes
out.Ne = Ne;          % Number of exposed nodes
out.Ni = Ni;          % Number of infected nodes
out.Nrec = Nrec;      % Number of recovered
out.NreExp = NreExp;  % number of nodes that are reexposed
out.NreInf = NreInf;  % number of nodes that are reinfected
out.NreRec = NreRec;  % number of nodes that are re-recovered
out.X = XX;           % Curent State 


end

