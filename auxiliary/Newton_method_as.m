%newton method with hessian and gradient of log posterior using asymmetric noise
function [minimizer] =newton_method_as(theta,data,mean,k,K,alpha,beta,n,m,lam,theta_mean,N)
for i =1:N
    theta = theta - (Asymmetric_hess_log(theta,data,mean,k,K,alpha,beta,n,m,lam)\eye((n+m)*n))*Asymmetric_grad_log(theta,data,mean,k,K,alpha,beta,n,m,lam,theta_mean);
end
minimizer = theta;
end