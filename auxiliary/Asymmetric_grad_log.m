%gradient of log posterior with asymmetric noise
function[grad_sum] = grad_log_as(theta,data,mean,k,K,alpha,beta,n,m,lam,theta_mean)

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

    %to satifsy the noise condition that expectation is zero,
    %translate the original noise by mean of it
    w=x_prime-Theta'*z;
    w_ = w+mean; 

    %calculate the gradient of log likelihood
    grad_log_theta = -Asymmetric_grad_log_w(w_,k,K,alpha,beta,n,1)*z;
    for j=2:n
        grad_log_theta = cat(2, grad_log_theta',(-Asymmetric_grad_log_w(w_,k,K,alpha,beta,n,j)*z)')';
    end

    grad_sum = grad_sum+grad_log_theta;
    
end