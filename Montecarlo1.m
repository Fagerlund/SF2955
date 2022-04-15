clear all; close all; clc
%Part 1
rng(19)
sigma=0.5; deltaT=0.5; alpha = 0.6;
P=1/20*[[16,1,1,1,1];[1,16,1,1,1];[1,1,16,1,1];[1,1,1,16,1];[1,1,1,1,16]];
zn=[[0,0];[3.5,0];[0,3.5];[0,-3.5];[-3.5,0]];
theta=[[1,deltaT,deltaT^2/2];[0,1,deltaT];[0,0,alpha]];
phiz=[[deltaT^2/2];[deltaT];[0]];
phiw=[[deltaT^2/2];[deltaT];[1]];
N=6;
zero=zeros(N,1);

N=6;
sigmamatrix=diag([500,5,5,200,5,5]);
%sigmamatrix=500*diag(ones(N,1),0) + 5*diag(ones(N-1,1),1) + 5*diag(ones(N-1,1),-1)+ 5*diag(ones(N-2,1),-2)+ 5*diag(ones(N-2,1),2)+ 200*diag(ones(N-3,1),3)+ 200*diag(ones(N-3,1),-3)+ 5*diag(ones(N-4,1),-4)+ 5*diag(ones(N-4,1),4)+ 5*diag(ones(N-5,1),5)+ 5*diag(ones(N-5,1),-5);

theta=[[theta,zeros(3,3)];[zeros(3,3),theta]];
phiz=[[phiz,zeros(3,1)];[zeros(3,1),phiz]];
phiw=[[phiw,zeros(3,1)];[zeros(3,1),phiw]];

X0=mvnrnd(zero,sigmamatrix);
Command=round(rand(1)*4)+1;
Z0=zn(Command,:);


m=3000; % Timesteps
States=zeros(m,2);
Xi=transpose(X0);
Zi=transpose(Z0);
Wi=transpose(mvnrnd([0,0],sigma^2*eye(2)));
for i=1:m
    
    Xi=theta*(Xi)+phiz*(Zi)+phiw*(Wi); %Calculate each new state
    
    %Need to update Zi,
    Update=(rand(1));
    sum=0;
    for j=1:5
        
        sum=P(Command,j)+sum;
        if Update<sum
            Command=j;
            break
        end
        
    end
    Zi=transpose(zn(Command,:));
    States(i,1)=Xi(1);
    States(i,2)=Xi(4);
    Wi=transpose(mvnrnd([0,0],sigma^2*eye(2)));
    
    
end
plot(States(:,1),States(:,2))
title('Randomized path, seed 19, m=3000')
xlabel('X1 location') 
ylabel('X2 location') 


%% Part 2

