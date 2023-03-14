%newton method with hessian and gradient of log posterior using gaussian mixture noise
function [minimizer] =newton_method_gm(theta,data,gm_a,r,n,m,lam,theta_mean,N)
for i =1:N
    theta = theta - (Gaussian_mixture_hess_log(theta,data,gm_a,r,n,m,lam)\eye((n+m)*n))*Gaussian_mixture_grad_log(theta,data,gm_a,r,n,m,lam,theta_mean);
end
minimizer = theta;
end