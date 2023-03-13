%hessian of log p(w)
function [a] = hess(w,k,K,alpha,beta,n)

hess_log_w = -k*eye(n);

if (alpha < w(n) <beta)
    hess_log_w(n,n) = -((K-k)/(beta-alpha)*(w(n))+(k-(K-k)*alpha/(beta-alpha)));

elseif (w(n)>=beta)
    hess_log_w(n,n) = -K;

end

a = hess_log_w;
