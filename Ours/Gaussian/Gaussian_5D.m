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

%Initial prior
lam = 1; %lambda
Theta_mean = zeros(n+m,n); %mean of the prior 

%Two alternate control
K1 = -eye(n);
K2 = randn(n);

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
    
    %calculate the mean of posterior(closed form)
    precision = eye(n+m); %inverse of covariance
    Theta_mean_cal = lam*Theta_mean;
    
    for k = 1:length(data)
        x = data{k}{1};
        u = data{k}{2};
        x_prime = data{k}{3};
        z = cat(2, x',u')';
        zeta = z* z';
        for i =1:n
            Theta_mean_cal(:,i) = Theta_mean_cal(:,i)+  x_prime(i)*z;
    
        end
        precision = precision + zeta;
    end
    sigma = 1*(precision)\eye(n+m); %covariance

    for l =1:n
        Theta_mean(:,l) = sigma*Theta_mean_cal(:,l);
    end
    

    disp(Theta_mean-Theta_true)
    error(t/interval,1) = norm(Theta_mean-Theta_true)/norm(Theta_true);
end
figure
hold on
x = linspace(interval,T,T/interval);
plot(x,error(1:T/interval),'r-.')
leg = legend('Gaussian');
set(leg,'Fontsize',10)

xlabel('Horizon','Fontsize',16)
ylabel('$|\theta-\theta_*|/|\theta_*|$','Interpreter','latex')
FileName = [datestr(now,'mmmm dd, yyyy HH:MM:SS.FFF AM')];
writematrix(error,['Ours_gaussian_5D',FileName '.csv']);
