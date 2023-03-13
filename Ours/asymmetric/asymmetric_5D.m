n=5;
m=5;

T=300;
N=200;

%{
T=50;
N=100;
%}

interval = 10;
error = zeros(T/interval,1);

%Stabilizable system parameters
A = [[0.5 0 0.6 0.3 0.1];[0.2 0 0.1 0.4 0.2];[0 0.2 0 0.5 0.3];[0.6 0.1 0.2 0 0.3];[0.1 0.4 0.3 0.5 0.1]];
B = [[0.4 0 0.5 0 0.4];[0 0.1 0.2 0.3 0.1];[0.3 0 0 0.2 0.1];[0.1 0.5 0.1 0.2 0.4];[0.4 0 0.5 0.3 0.1]];

%Unstabilizable system parameters
%{
A=[[1 0 0 0 0];[1 0 1 0 -1];[0 0 -1 0 1];[1 1 1 0 0];[1 0 -1 0 0]];
B=[[0 0 0 0 0];[1 0 0 0 0];[0 0 0 0 0];[1 0 1 0 0];[1 0 0 0 0]];
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

%Asymmetric noise setting
k=1;
K=10;
alpha = 0;
beta = 40;

%To satifsy the noise condition that expectation is zero,
%calculate the mean of original noise and substract it.
noise = readmatrix('asymmetric_noise_1D.csv');
mean = [0 0 0 0 -0.17873]';

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
            
            rn=randi(100000);
            as_noise = [randn(n-1,1)',noise(rn)]'-mean;
            x_prime = A*x + B*u + as_noise;

            data{end+1} = {x;u;x_prime};
            x=x_prime;
        end
        state_data(i,:) = x';
    end

    %argmax of U by newton method 
    
    theta = newton_method(theta,data,mean,k,K,alpha,beta,n,m,lam,theta_mean,100);    

    disp(theta-theta_true)
    error(t/interval,1) = norm(theta-theta_true)/norm(theta_true);
end
figure
hold on
x = linspace(interval,T,T/interval);
plot(x,distance(1:T/interval),'r-.')
leg = legend('asymmetric');
set(leg,'Fontsize',10)

xlabel('Horizon','Fontsize',16)
ylabel('$|\theta-\theta_*|$','Interpreter','latex')
FileName = [datestr(now,'mmmm dd, yyyy HH:MM:SS.FFF AM')];
writematrix(distance,['ours_asymmetric_5D',FileName '.csv']);

%newton method with hessian and gradient of log posterior using asymmetric noise
function [minimizer] =newton_method(theta,data,mean,k,K,alpha,beta,n,m,lam,theta_mean,N)
for i =1:N
    theta = theta - (asymmetric_hess_log(theta,data,mean,k,K,alpha,beta,n,m,lam)\eye((n+m)*n))*asymmetric_grad_log(theta,data,mean,k,K,alpha,beta,n,m,lam,theta_mean);
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
