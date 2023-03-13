function[hess_sum] = hess_log_as(theta,data,mean,k,K,alpha,beta,n,m,lam)

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

    %to satifsy the noise condition that expectation is zero,
    %translate the original noise by mean of it
    w=x_prime-Theta'*z;
    w_ = w+mean; 

    %calculate the hessian of log likelihood
    hess_log_phi = kron(hess_log_w(w_,k,K,alpha,beta,n),z*z');

    hess_sum = hess_sum+hess_log_phi;
    
end