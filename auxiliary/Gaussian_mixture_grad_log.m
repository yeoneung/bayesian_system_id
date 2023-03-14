%gradient of log posterior with gaussian mixture noise
function[grad_sum] = grad_log_gm(theta,data,gm_a,r,n,m,lam,theta_mean)

grad_sum = -lam*(theta-theta_mean); %intialize the gradient of log posterior as gradient of log prior

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
    
    %calculate the gradient of log likelihood
    grad_log_theta = 1/r*(w(1)-gm_a(1)+2*gm_a(1)/(1+exp(1)^(2/r*w'*gm_a)))*z;
    for j=2:n
        grad_log_theta = cat(2, grad_log_theta',(1/r*(w(j)-gm_a(j)+2*gm_a(j)/(1+exp(1)^(2/r*w'*gm_a)))*z)')';
    end

    grad_sum = grad_sum+grad_log_theta;
    
end