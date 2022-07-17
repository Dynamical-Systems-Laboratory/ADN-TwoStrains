
function xxi = activityPotential(alfa,xmin,xmax,N)

for i=1:N
    
    r=rand();
    xxi(i)=power((power(xmax,alfa+1)-power(xmin,alfa+1))*r+power(xmin,alfa+1),(1/(alfa+1)));
end
