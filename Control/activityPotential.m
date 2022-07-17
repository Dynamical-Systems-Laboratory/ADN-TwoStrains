function xxi = activityPotential(alfa,xmin,xmax,N)

% F = @(x)x.^-gamma;
% 
% x = 1 + (5-1)*rand(N,1); % uniform random numbers on [e,1]
% xi = F(x);
% epsilon = 1e-3;
% 
% xi(xi<epsilon) = epsilon;
% 
% a = eta*xi;  % Activity

%%Epidemic threshold


for i=1:N
    
    r=rand();
    xxi(i)=power((power(xmax,alfa+1)-power(xmin,alfa+1))*r+power(xmin,alfa+1),(1/(alfa+1)));
end
