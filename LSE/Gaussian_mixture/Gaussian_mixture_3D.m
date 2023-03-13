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

% Gaussian mixture setting
gm_a = [1/4, 1/4, 1/4]';
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

            %randomized action 
            u = mvnrnd(zeros(n,1),1*eye(n))';

            %Gaussian noise
            w = random(gm)';           

            x_prime = A*x + B*u + w;
                       
            if j ==interval              
                data{end+1} = {x;u;x_prime;w};
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
leg = legend('Gaussian_mixture');
set(leg,'Fontsize',10)

xlabel('Horizon','Fontsize',16)
ylabel('$|\theta-\theta_*|/|\theta_*|$','Interpreter','latex')
FileName = [datestr(now,'mmmm dd, yyyy HH:MM:SS.FFF AM')];
writematrix(error,['LSE_gaussian_mixture_3D',FileName '.csv']);
