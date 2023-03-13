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

%Asymmetric noise setting
M=1;
L=10;
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

            %randomized action 
            u = mvnrnd(zeros(n,1),1*eye(n))';
            
            %asymmetric noise
            rn=randi(100000);
            as_noise = [randn(4,1)',noise(rn)]'-mean;

            x_prime = A*x + B*u + as_noise;
                       
            if j ==interval              
                data{end+1} = {x;u;x_prime;as_noise};
            end            
            x = x_prime;
        end
        state_data(i,:) = x';
    end
    
    %argmax of U  
    precision = 0;
    Theta_mean_cal = zeros(n+m,n);
    Theta_mean = zeros(n+m,n);
    
    for k = 1:length(data)
        x = data{k}{1};
        u = data{k}{2};
        x_prime = data{k}{3};
        w = data{k}{4};
        z = cat(2, x',u')';
        zeta = z* z';
        for i =1:n
            Theta_mean_cal(:,i) = Theta_mean_cal(:,i)+  w(i)*z; 
        end
        precision = precision + zeta;
    end
    sigma = 1*(precision)\eye(n+m);

    for l =1:n
        Theta_mean(:,l) = sigma*Theta_mean_cal(:,l);
    end
    
    %initialize data set
    data = {};

    disp(Theta_mean)
    error(t/interval,1) = norm(Theta_mean)/norm(Theta_true);
end
figure
hold on
x = linspace(interval,T,T/interval);
plot(x,error(1:T/interval),'r-.')
leg = legend('asymmetric');
set(leg,'Fontsize',10)

xlabel('Horizon','Fontsize',16)
ylabel('$|\theta-\theta_*|/|\theta_*|$','Interpreter','latex')
FileName = [datestr(now,'mmmm dd, yyyy HH:MM:SS.FFF AM')];
writematrix(error,['LSE_asymmetric_5D',FileName '.csv']);
