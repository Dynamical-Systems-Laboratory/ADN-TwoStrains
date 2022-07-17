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

function [Lambda,Lambda2] = stability_switchedSys(Parameters,Par1,Par2)

T = Parameters.Wtime;
D = Parameters.WDuration;
p = Parameters.p;
dt = Parameters.dt;
m = Parameters.m;
a = Parameters.a;

am = mean(a);
am2 = mean(a.^2);

for p1=1:length(Par1)

    Parameters.al11 = Par1(p1);
    Parameters.al22 = Par1(p1);
    

    for p2=1:length(Par2)


        Parameters.lambda1 = Par2(p2);
        Parameters.lambda2 = Par2(p2);



        mu = Parameters.mu1;
        sigma = Parameters.sigma_1;

        % Switching times
        t1 = D*T*dt;
        t2 = (1-D)*T*dt;

        %% stability for epidemic equilibrium
        c1 = m*Parameters.lambda1;
        J0 = [ -mu,       0,     sigma,        0;
            0,       -mu,      0,       sigma;
            c1*am,       c1,   -sigma,       0;
            c1*am2,      c1*am,    0,      -sigma ];

        c1 = m*Parameters.lambda1*(1-p/100);

        Jp = [ -mu,       0,    sigma,        0;
            0,       -mu,      0,       sigma;
            c1*am,      c1,  -sigma,       0;
            c1*am2,   c1*am,    0,       -sigma];

        Phi = expm(J0*t2)*expm(Jp*t1);

        [R,~] = logm(Phi);
        R = (1/(t1+t2)) * R;
        Eig_m= eig(R);

        Lambda(p1,p2) = max(real(Eig_m));

        %% stability for endemic equilibrium
        mu = Parameters.mu1;
        sigma = Parameters.sigma_1;
        c1 = m*Parameters.al11*Parameters.lambda1;

        J0 = [ -mu,       0,    sigma,        0;
            0,       -mu,      0,       sigma;
            c1*am,      c1,  -sigma,       0;
            c1*am2,   c1*am,    0,      -sigma ];

        c1 = m*Parameters.al11*Parameters.lambda1*(1-p/100);

        Jp = [ -mu,       0,    sigma,        0;
                0,       -mu,      0,       sigma;
              c1*am,      c1,  -sigma,       0;
              c1*am2,   c1*am,    0,       -sigma];

        Phi = expm(J0*t2)*expm(Jp*t1);

        [R,~] = logm(Phi);
        R = 1/(t1+t2)*R;
        Eig_m= eig(R);
        Lambda2(p1,p2) = max(real(Eig_m));

    end

end
end
