n=3;
m=3;

T=300;
N=200;
interval = 10;
error = zeros(T/interval,1);

%Stabilizable system parameters
A = [[0.5 0 0.6];[0.2 0 0.1];[0 0.5 0.3]];
B = [[0.4 0 0.4];[0 0.3 0.1];[0.3 0.2 0.1]];

%Unstabilizable system parameters
%{
A=[[0 -1 1];[1 0 1];[0 0 1]];
B=[[1 0 0];[0 0 1];[0 0 0]];
%}

Theta_true = [A,B]';
theta_true = tr_Theta_to_theta(Theta_true,n,m); % vectorized true system parameter

%Initial prior
lam = 1; %lambda
Theta_mean = zeros(n+m,n); %mean of the prior 
theta_mean = tr_Theta_to_theta(Theta_mean,n,m); 
theta = theta_mean; % initial vectorized system parameter

%Two alternate control
K1 = -eye(n);
K2 = randn(n);

% Gaussian mixture setting
gm_a = [1/4,1/4,1/4]';
mu = [gm_a';-gm_a'];
r=1;
cov = [r,r,r];
gm = gmdistribution(mu,cov);

%we observe errors every interval
interval=10;
data = {};
state_data = zeros(N,n);

%main algorithm
for t = interval:interval:T
    disp(t)
    for i = 1:N
        x=(state_data(i,:))';
        for j = 1:interval 
             
            u = K1*x;
    
            if rem(j,2)==0
                
                u = K2*x;
            end
            
            x_prime = A*x + B*u +mvnrnd(zeros(n,1),1*eye(n))';
            data{end+1} = {x;u;x_prime};
            x=x_prime;
        end
        state_data(i,:) = x';
    end

    %argmax of U by newton method 
    
    theta = newton_method(theta,data,gm_a,r,n,m,lam,theta_mean,100);    

    disp(theta-theta_true)
    error(t/interval,1) = norm(theta-theta_true)/norm(theta_true);
end
figure
hold on
x = linspace(interval,T,T/interval);
plot(x,distance(1:T/interval),'r-.')
leg = legend('Gaussian mixture');
set(leg,'Fontsize',10)

xlabel('Horizon','Fontsize',16)
ylabel('$|\theta-\theta_*|$','Interpreter','latex')
FileName = [datestr(now,'mmmm dd, yyyy HH:MM:SS.FFF AM')];
writematrix(distance,['ours_gaussian_mixture_3D',FileName '.csv']);


%newton method with hessian and gradient of log posterior using gaussian mixture noise
function [minimizer] =newton_method(theta,data,gm_a,r,n,m,lam,theta_mean,N)
for i =1:N
    theta = theta - (gaussian_mixture_hess_log(theta,data,gm_a,r,n,m,lam)\eye((n+m)*n))*gaussian_mixture_grad_log(theta,data,gm_a,r,n,m,lam,theta_mean);
end
minimizer = theta;
end

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

