%stabilizable system parameters 
A_stab_3D = [[0.5 0 0.6];[0.2 0 0.1];[0 0.5 0.3]];
B_stab_3D = [[0.4 0 0.4];[0 0.3 0.1];[0.3 0.2 0.1]];

A_stab_5D = [[0.5 0 0.6 0.3 0.1];[0.2 0 0.1 0.4 0.2];[0 0.2 0 0.5 0.3];[0.6 0.1 0.2 0 0.3];[0.1 0.4 0.3 0.5 0.1]];
B_stab_5D = [[0.4 0 0.5 0 0.4];[0 0.1 0.2 0.3 0.1];[0.3 0 0 0.2 0.1];[0.1 0.5 0.1 0.2 0.4];[0.4 0 0.5 0.3 0.1]];

%unstabilizable system parameters
A_unstab_3D = [[0 -1 1];[1 0 1];[0 0 1]];
B_unstab_3D = [[1 0 0];[0 0 1];[0 0 0]];

A_unstab_5D = [[1 0 0 0 0];[1 0 1 0 -1];[0 0 -1 0 1];[1 1 1 0 0];[1 0 -1 0 0]];
B_unstab_5D = [[0 0 0 0 0];[1 0 0 0 0];[0 0 0 0 0];[1 0 1 0 0];[1 0 0 0 0]];

%Initial prior
lam = 1; %lambda
Theta_mean_3D = zeros(6,3); %mean of the prior in 3D
Theta_mean_5D = zeros(10,5); %mean of the prior in 5D

%stabilizable system paratmer with Gaussian noise
%n=3, m=3, T=300, N=200, interval=10
Our_Gaussian(A_stab_3D,B_stab_3D,lam,Theta_mean_3D,3,3,300,200,10)

%n=5, m=5, T=300, N=200, interval=10
Our_Gaussian(A_stab_5D,B_stab_5D,lam,Theta_mean_5D,5,5,300,200,10)

%n=3, m=3, T=300, N=200, interval=10
LSE_Gaussian(A_stab_3D,B_stab_3D,3,3,300,200,10)

%n=5, m=5, T=300, N=200, interval=10
LSE_Gaussian(A_stab_5D,B_stab_5D,5,5,300,200,10)


%stabilizable system paratmer with Gaussian mixture noise
%n=3, m=3, T=300, N=200, interval=10
Our_Gaussian_mixture(A_stab_3D,B_stab_3D,lam,Theta_mean_3D,3,3,300,200,10)

%n=5, m=5, T=300, N=200, interval=10
Our_Gaussian_mixture(A_stab_5D,B_stab_5D,lam,Theta_mean_5D,5,5,300,200,10)

%n=3, m=3, T=300, N=200, interval=10
LSE_Gaussian_mixture(A_stab_3D,B_stab_3D,3,3,300,200,10)

%n=5, m=5, T=300, N=200, interval=10
LSE_Gaussian_mixture(A_stab_5D,B_stab_5D,5,5,300,200,10)


%stabilizable system paratmer with asymmetric noise
%n=3, m=3, T=300, N=200, interval=10
Our_asymmetric(A_stab_3D,B_stab_3D,lam,Theta_mean_3D,3,3,300,200,10)

%n=5, m=5, T=300, N=200, interval=10
Our_asymmetric(A_stab_5D,B_stab_5D,lam,Theta_mean_5D,5,5,300,200,10)

%n=3, m=3, T=300, N=200, interval=10
LSE_asymmetric(A_stab_3D,B_stab_3D,3,3,300,200,10)

%n=5, m=5, T=300, N=200, interval=10
LSE_asymmetric(A_stab_5D,B_stab_5D,5,5,300,200,10)


%unstabilizable system paratmer with Gaussian noise
%n=3, m=3, T=50, N=100, interval=10
Our_Gaussian(A_unstab_3D,B_unstab_3D,lam,Theta_mean_3D,3,3,50,100,10)

%n=5, m=5, T=50, N=100, interval=10
Our_Gaussian(A_unstab_5D,B_unstab_5D,lam,Theta_mean_5D,5,5,50,100,10)

%n=3, m=3, T=50, N=100, interval=10
LSE_Gaussian(A_unstab_3D,B_unstab_3D,3,3,50,100,10)

%n=5, m=5, T=50, N=100, interval=10
LSE_Gaussian(A_unstab_5D,B_unstab_5D,5,5,50,100,10)


%unstabilizable system paratmer with Gaussian mixture noise
%n=3, m=3, T=50, N=100, interval=10
Our_Gaussian_mixture(A_unstab_3D,B_unstab_3D,lam,Theta_mean_3D,3,3,50,100,10)

%n=5, m=5, T=50, N=100, interval=10
Our_Gaussian_mixture(A_unstab_5D,B_unstab_5D,lam,Theta_mean_5D,5,5,50,100,10)

%n=3, m=3, T=50, N=100, interval=10
LSE_Gaussian_mixture(A_unstab_3D,B_unstab_3D,3,3,50,100,10)

%n=5, m=5, T=50, N=100, interval=10
LSE_Gaussian_mixture(A_unstab_5D,B_unstab_5D,5,5,50,100,10)


%unstabilizable system paratmer with asymmetric noise
%n=3, m=3, T=50, N=100, interval=10
Our_asymmetric(A_unstab_3D,B_unstab_3D,lam,Theta_mean_3D,3,3,50,100,10)

%n=5, m=5, T=50, N=100, interval=10
Our_asymmetric(A_unstab_5D,B_unstab_5D,lam,Theta_mean_5D,5,5,50,100,10)

%n=3, m=3, T=50, N=100, interval=10
LSE_asymmetric(A_unstab_3D,B_unstab_3D,3,3,50,100,10)

%n=5, m=5, T=50, N=100, interval=10
LSE_asymmetric(A_unstab_5D,B_unstab_5D,5,5,50,100,10)

