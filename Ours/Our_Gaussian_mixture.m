function[] = ours_gaussian_mixture(A,B,lam,Theta_mean,n,m,T,N,interval)
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

% Gaussian mixture setting
if n==3
    gm_a = [1/4,1/4,1/4]';
    mu = [gm_a';-gm_a'];
    r=1;
    cov = [r,r,r];
else
    gm_a = [1/4,1/4,1/4,1/4,1/4]';
    mu = [gm_a';-gm_a'];
    r=1;
    cov = [r,r,r,r,r];
end
gm = gmdistribution(mu,cov);

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
            
            x_prime = A*x + B*u +random(gm)';
            data{end+1} = {x;u;x_prime};
            x=x_prime;
        end
        state_data(i,:) = x';
    end

    %argmax of U by newton method 
    
    theta = Newton_method_gm(theta,data,gm_a,r,n,m,lam,theta_mean,100);    

    disp(theta-theta_true)
    error(t/interval,1) = norm(theta-theta_true)/norm(theta_true);
end
figure
hold on
x = linspace(interval,T,T/interval);
plot(x,error(1:T/interval),'r-.')
leg = legend('Gaussian mixture');
set(leg,'Fontsize',10)

xlabel('Horizon','Fontsize',16)
ylabel('$|\theta-\theta_*|$','Interpreter','latex')
save_time = [datestr(now,'mmmm dd, yyyy HH:MM:SS.FFF AM')];
FileName = "Our_Gaussian_mixture-"+string(n)+"D_"+save_time+".csv";
writematrix(error,'FileName');
end
