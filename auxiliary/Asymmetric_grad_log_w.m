%i-th component of gradient log p(w)
function [a] = grad(w,k,K,alpha,beta,n,i)
C_1 = (K-k)/2*alpha^2/(beta-alpha);
D_1 = -(K-k)/2*(alpha+beta);

grad_log_w = -k*w;

if (alpha < w(n) <beta)
    grad_log_w(n) = -((K-k)/(2*(beta-alpha))*(w(n))^2+(k-(K-k)*alpha/(beta-alpha))*w(n)+C_1);

elseif (w(n)>=beta)
    grad_log_w(n) = -(K*w(n)+D_1);

end

a = grad_log_w(i);
