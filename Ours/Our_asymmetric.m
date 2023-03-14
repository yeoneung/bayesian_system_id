function[] = ours_asymmetric(A,B,lam,Theta_mean,n,m,T,N,interval)
%n = state dimension
%m = control dimension
%T = Time horizon
%N = rollout
%we observe error between true system parmaeter and estimator every interval
error = zeros(T/interval,1);

Theta_true = [A,B]';
theta_true = Tr_Theta_to_theta(Theta_true,n,m); % vectorized true system parameter

theta_mean = Tr_Theta_to_theta(Theta_mean,n,m); %vectorized mean of the prior
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
if n==3
    mean = [0 0 -0.17873]';
else
    mean = [0 0 0 0 -0.17873]';
end

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
    
    theta = Newton_method_as(theta,data,mean,k,K,alpha,beta,n,m,lam,theta_mean,100);    

    disp(theta-theta_true)
    error(t/interval,1) = norm(theta-theta_true)/norm(theta_true);
end
figure
hold on
x = linspace(interval,T,T/interval);
plot(x,error(1:T/interval),'r-.')
leg = legend('asymmetric');
set(leg,'Fontsize',10)

xlabel('Horizon','Fontsize',16)
ylabel('$|\theta-\theta_*|$','Interpreter','latex')
save_time = [datestr(now,'mmmm dd, yyyy HH:MM:SS.FFF AM')];
FileName = "Our_asymmetric-"+string(n)+"D_"+save_time+".csv";
writematrix(error,FileName);
end
