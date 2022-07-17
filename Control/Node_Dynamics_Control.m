function [Xo] = Node_Dynamics_Control(Act,IDcontrol,X,Xn,Parameters)

m = Parameters.m;
a11 = Parameters.al11;
a12 = Parameters.al12;
a21 = Parameters.al21;
a22 = Parameters.al22;
lambda1 = Parameters.lambda1;
lambda2 = Parameters.lambda2;
sigma1 = Parameters.sigma_1;
sigma2 = Parameters.sigma_2;
mu1 = Parameters.mu1;
mu2 = Parameters.mu2;

dt = Parameters.dt;
N = Parameters.N;


Xo = Xn;

xrem = setdiff(1:N, IDcontrol');
if isempty(IDcontrol)==1
    nR = 0;
else
    nR = length(IDcontrol);
end


for i=1:N

    if X(i) == 1  % if node is exposed with variant 1 it can get infected with I1
        if rand<sigma1*dt
            Xo(i) = 2;
        end
        %
    elseif X(i) == 1.5  % if node is exposed with variant 2 it can get infected with I2
        if rand<sigma2*dt
            Xo(i) = 2.5;
        end
        %
    elseif X(i) == 2    % If node is infected with I1 then it can die or recover
        if rand<mu1*dt
            Xo(i) = 3;  % Node recovers
        end
        %
        %
    elseif X(i) == 2.5    % If node is infected with I2 then it can die or recover
        if rand<mu2*dt
            Xo(i) = 3.5;
        end
        %
        %
        %% reinfection dynamics
    elseif X(i) == 4      % E12: if node is exposed to variant 2 given it was infected from 1, it can get infected with I2
        if rand<sigma2*dt
            Xo(i) = 5;    % I12: infected with 2 given it was infected with 1 already
        end
        %
    elseif X(i) == 4.5    % E21 if node is exposed to variant 1 given it was infected from 2, it can get infected with I1
        if rand<sigma1*dt
            Xo(i) = 5.5;  % I21: infected with 1 given it was infected with 2 already
        end
        %
    elseif X(i) == 5      % If node is infected with I12 then it can die or recover
        if rand<mu2*dt
            Xo(i) = 6;    % Node recovers
        end
        %
        %
    elseif X(i) == 5.5    % If node is infected with I21 then it can die or recover
        if rand<mu1*dt
            Xo(i) = 6;
        end
    end



    if Act(i) == 1 %% Node becomes active
        ranID = 1+(N-nR-1)*rand(1,m);
        links = xrem(round(ranID)); % random indices neglecting selfloops

        for j=1:m
            l=links(j);

            % Transitions from the Suceptible state
            %% Transitions from S to E1
            % by contact with I1
            %  Xi = I1 and Xj = S
            if X(i) == 2 && X(l) == 0
                if rand<lambda1
                    Xo(l) = 1;                  %Xj = E1
                end
                %  Xi = S and Xj = I1
            elseif X(i) == 0 && X(l) == 2
                if rand<lambda1
                    Xo(i) = 1;                  %Xi = E1
                end
                %
                % by contact with I21
                %  Xi = I21 and Xj = S
            elseif X(i) == 5.5 && X(l) == 0
                if rand<lambda1
                    Xo(l) = 1;                  %Xj = E1
                end
                %  Xi = S and Xj = I21
            elseif X(i) == 0 && X(l) == 5.5
                if rand<lambda1
                    Xo(i) = 1;                  %Xi = E1
                end
                %
                %% Transitions from S to E2
                % by contact with I2
                %  Xi = I2 and Xj = S
            elseif X(i) == 2.5 && X(l) == 0
                if rand<lambda2
                    Xo(l) = 1.5;                  %Xj = E2
                end
                %  Xi = S and Xj = I2
            elseif X(i) == 0 && X(l) == 2.5
                if rand<lambda2
                    Xo(i) = 1.5;                  %Xi = E2
                end
                %
                % by contact with I12
                %  Xi = I12 and Xj = S
            elseif X(i) == 5 && X(l) == 0
                if rand<lambda2
                    Xo(l) = 1.5;                   %Xj = E2
                end
                %  Xi = S and Xj = I21
            elseif X(i) == 0 && X(l) == 5
                if rand<lambda2
                    Xo(i) = 1.5;                  %Xi = E2
                end
                %
                % Transitions from the Recovered states R1 and R2
                %% Transitions from R1 to E1
                % by contact with I1
                %  Xi = I1 and Xj = R1: Reinfection
            elseif X(i) == 2 && X(l) == 3
                if rand<a11*lambda1
                    Xo(l) = 1;                  %Xj = E1
                end
                %  Xi = R1 and Xj = I1: Reinfection
            elseif X(i) == 3 && X(l) == 2
                if rand<a11*lambda1
                    Xo(i) = 1;                  %Xj = E1
                end
                % by contact with I21
                % Xi = I21 and Xj = R1: Reinfection
            elseif X(i) == 5.5 && X(l) == 3
                if rand<a11*lambda1
                    Xo(l) = 1;                  %Xj = E1
                end
                %  Xi = R1 and Xj = I21: Reinfection
            elseif X(i) == 3 && X(l) == 5.5
                if rand<a11*lambda1
                    Xo(i) = 1;                  %Xj = E1
                end
                %
                %% Transitions from R1 to E12
                % by contact with I2
                %  Xi = I2 and Xj = R1: Reinfection
            elseif X(i) == 2.5 && X(l) == 3
                if rand<a12*lambda2
                    Xo(l) = 4;                  %Xj = E12
                end
                %  Xi = R1 and Xj = I1: Reinfection
            elseif X(i) == 3 && X(l) == 2.5
                if rand<a12*lambda2
                    Xo(i) = 4;                  %Xj = E12
                end
                % by contact with I12
                %  Xi = I12 and Xj = R1: Reinfection
            elseif X(i) == 5 && X(l) == 3
                if rand<a12*lambda2
                    Xo(l) = 4;                  %Xj = E12
                end
                %  Xi = R1 and Xj = I1: Reinfection
            elseif X(i) == 3 && X(l) == 5
                if rand<a12*lambda2
                    Xo(i) = 4;                  %Xj = E12
                end
                %% Transitions from R2 to E2
                %  by contact with I2
                %  Xi = I2 and Xj = R2: Reinfection
            elseif X(i) == 2.5 && X(l) == 3.5
                if rand<a22*lambda2
                    Xo(l) = 1.5;                  %Xj = E2
                end
                % Xi = R2 and Xj = I2: Reinfection
            elseif X(i) == 3.5 && X(l) == 2.5
                if rand<a22*lambda2
                    Xo(i) = 1.5;                  %Xj = E2
                end
                % by contact with I21
                % Xi = I21 and Xj = R2: Reinfection
            elseif X(i) == 5 && X(l) == 3.5
                if rand<a22*lambda2
                    Xo(l) = 1.5;                  %Xj = E2
                end
                %  Xi = R1 and Xj = I21: Reinfection
            elseif X(i) == 3.5 && X(l) == 5
                if rand<a22*lambda2
                    Xo(i) = 1.5;                  %Xj = E2
                end
                %
                %% Transitions from R2 to E21
                %  by contact with I1
                %  Xi = I1 and Xj = R2: Reinfection
            elseif X(i) == 2 && X(l) == 3.5
                if rand<a21*lambda1
                    Xo(l) = 4.5;                  %Xj = E21
                end
                %  Xi = R2 and Xj = I1: Reinfection
            elseif X(i) == 3.5 && X(l) == 2
                if rand<a21*lambda1
                    Xo(i) = 4.5;                  %Xj = E21
                end
                %  by contact with I21
                %  Xi = I21 and Xj = R2: Reinfection
            elseif X(i) == 5.5 && X(l) == 3.5
                if rand<a21*lambda1
                    Xo(l) = 4.5;                  %Xj = E21
                end
                %  Xi = R2 and Xj = I21: Reinfection
            elseif X(i) == 3.5 && X(l) == 5.5
                if rand<a21*lambda1
                    Xo(i) = 4.5;                  %Xj = E21
                end
                %% Transitions from R to E12
                %  by contact with I2
                %  Xi = I2 and Xj = R: Reinfection
            elseif X(i) == 2.5 && X(l) == 6
                if rand<a22*lambda2
                    Xo(l) = 4;                  %Xj = E12
                end
            elseif X(i) == 6 && X(l) == 2.5
                if rand<a22*lambda2
                    Xo(i) = 4;                  %Xj = E12
                end
                %  by contact with I12
                %  Xi = I12 and Xj = R: Reinfection
            elseif X(i) == 5 && X(l) == 6
                if rand<a22*lambda2
                    Xo(l) = 4;                  %Xj = E12
                end
            elseif X(i) == 6 && X(l) == 5
                if rand<a22*lambda2
                    Xo(i) = 4;                  %Xj = E12
                end
                %
                %
                %% Transitions from R to E21
                %  by contact with I1
                %  Xi = I1 and Xj = R: Reinfection
            elseif X(i) == 2 && X(l) == 6
                if rand<a11*lambda1
                    Xo(l) = 4.5;                  %Xj = E21
                end
            elseif X(i) == 6 && X(l) == 2
                if rand<a11*lambda1
                    Xo(i) = 4.5;                  %Xj = E21
                end
                %  by contact with I12
                %  Xi = I12 and Xj = R: Reinfection
            elseif X(i) == 5.5 && X(l) == 6
                if rand<a11*lambda1
                    Xo(l) = 4.5;                  %Xj = E21
                end
            elseif X(i) == 6 && X(l) == 5.5
                if rand<a11*lambda1
                    Xo(i) = 4.5;                  %Xj = E21
                end
            end
            %
            %
        end                                    % end for cycle for the m random links
    end
end
end
