%vectorize system parameter
function [theta]=tr_Theta_to_theta(Theta,n,m)
    theta_=[];
    for j=1:n
        for i=1:n+m
            theta_(end+1)=Theta(i,j);
        end
    end
    theta=theta_';
end
