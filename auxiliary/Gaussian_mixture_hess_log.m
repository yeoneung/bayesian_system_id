function[hess_sum] = hess_log_gm(theta,data,gm_a,r,n,m,lam)

hess_sum = -lam*eye((n+m)*n); %intialize the hessian of log posterior as hessian of log prior

for transition = 1:length(data)
    x = data{transition}{1};
    u = data{transition}{2};
    x_prime = data{transition}{3};
    
    z = cat(2, x', u')';
    
    %transfer vectorized system parameter to matrix form
    Theta=zeros((n+m),n);
    for i=1:n
        Theta(:,i) = theta((n+m)*(i-1)+1:(n+m)*i);
    end
    
    w=x_prime-Theta'*z;

    %calculate the hessian of log likelihood
    hess_log_theta = kron(-((1/r)*eye(n)-4*(1/(r^2))*(gm_a*gm_a')*(exp(2*w'*gm_a/r)/(1+exp(2*w'*gm_a/r))^2)),z*z');

    hess_sum = hess_sum+hess_log_theta;
    
end