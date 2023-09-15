function[] = LSE_asymmetric(A,B,n,m,T,N,interval)
%n = state dimension
%m = control dimension
%T = Time horizon
%N = rollout
%we observe error between true system parmaeter and estimator every interval
error = zeros(T/interval,1);

Theta_true = [A,B]';

%Asymmetric noise setting
M=1;
L=10;
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

            %randomized action 
            u = mvnrnd(zeros(n,1),1*eye(n))';
            
            %asymmetric noise
            rn=randi(100000);
            as_noise = [randn(n-1,1)',noise(rn)]'-mean;

            x_prime = A*x + B*u + as_noise;
                       
            if j ==interval              
                data{end+1} = {x;u;x_prime;as_noise};
            end            
            x = x_prime;
        end
        state_data(i,:) = x';
    end
    
    %least square estimator
    %to obtain error between true system parameter and LSE, we use equation (12) in Dean2020
    precision = 0;
    Error_cal = zeros(n+m,n);
    Error = zeros(n+m,n);
    
    for k = 1:length(data)
        x = data{k}{1};
        u = data{k}{2};
        x_prime = data{k}{3};
        w = data{k}{4};
        z = cat(2, x',u')';
        zeta = z* z';
        for i =1:n
            Error_cal(:,i) = Error_cal(:,i)+  w(i)*z;
    
        end
        precision = precision + zeta;
    end
    sigma = 1*(precision)\eye(n+m);

    for l =1:n
        Error(:,l) = sigma*Error_cal(:,l);
    end
    
    %initialize data set
    data = {};

    disp(Error)
    error(t/interval,1) = norm(Error)/norm(Theta_true);
end
figure
hold on
x = linspace(interval,T,T/interval);
plot(x,error(1:T/interval),'r-.')
leg = legend('asymmetric');
set(leg,'Fontsize',10)

xlabel('Horizon','Fontsize',16)
ylabel('$|\theta-\theta_*|/|\theta_*|$','Interpreter','latex')
save_time = [datestr(now,'mmmm dd, yyyy HH:MM:SS.FFF AM')];
FileName = "LSE_asymmetric-"+string(n)+"D_"+save_time+".csv";
writematrix(error,'FileName');
end
