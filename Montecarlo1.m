clear all; close all; clc
%Part 1
%rng(19)
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
Command=randi([1 5],1,1);
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
    summ=0;
    for j=1:5
        
        summ=P(Command,j)+summ;
        if Update<summ
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

%% Making of a nxN state matrix for use in part 3 
v=90; zeta=1.5; gamma=3; 

load("stations.mat");

load("RSSI-measurements.mat");

N=10000;
n=length(Y);
Command=randi([1 5],1,N);
Commands=zeros(N,n);
Commands(:,1)=Command;
for j=1:N
    
    for i =1:n
        Update=(rand(1));
        summ=0;
        for k=1:5

            summ=P(round(Commands(j,i)),k)+summ;
            if Update<=summ
                Command=k;
                break
            end

        end
        Commands(j,i+1)=Command;
         
    end
    
end

%% Part 3



tau=zeros(6,n);
Xi=transpose(mvnrnd(zero,sigmamatrix,N)); %Initilization

%p = @(x1,x2,y) prod(normpdf(Y,(6*v*ones(1,N)-10*gamma*log10(vecnorm(([x1;x2]-pos_vec(:,1).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,2).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,3).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,4).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,5).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,6).*ones(2,N))))),zeta),2);%1/(2*pi*zeta^2)^3*exp(-1/2*(sum(y)*ones(1,N)-(6*v*ones(1,N)-10*gamma*log10(vecnorm(([x1;x2]-pos_vec(:,1).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,2).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,3).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,4).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,5).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,6).*ones(2,N)))))));
w = (p(Xi(1,:),Xi(4,:),Y(:,1)));
tau(1,1) = sum(Xi(1,:).*w)/sum(w);
tau(2,1) = sum(Xi(4,:).*w)/sum(w);

Wi=transpose(mvnrnd([0,0],sigma^2*eye(2),N));


Zi=transpose(zn(Commands(:,1),:));
for k = 1:500 % main loop
    Xi=theta*(Xi)+phiz*(Zi)+phiw*(Wi); %Calculate each new state
    w =w.*(p(Xi(1,:),Xi(4,:),Y(:,k+1))); % weighting REMOVE/INCLUDE LOG
    
    %w=w/sum(w);
    Zi=transpose(zn(Commands(:,k+1),:));
    Wi=transpose(mvnrnd([0,0],sigma^2*eye(2),N));
    tau(1,k+1) =  sum(Xi(1,:).*w)/sum(w);
    tau(2,k+1) =  sum(Xi(4,:).*w)/sum(w);
end

plot(tau(1,:),tau(2,:))
title("SIS, N=10000")
hold on
plot(pos_vec(1,1),pos_vec(2,1),'d')
hold on
plot(pos_vec(1,2),pos_vec(2,2),'d')
hold on
plot(pos_vec(1,3),pos_vec(2,3),'d')
hold on
plot(pos_vec(1,4),pos_vec(2,4),'d')
hold on
plot(pos_vec(1,5),pos_vec(2,5),'d')
hold on
plot(pos_vec(1,6),pos_vec(2,6),'d')
%figure
%hist(w)
ESSM=1/(sum((w/sum(w)).^2))

%% Part 3 Instructed weight adjustment



tau=zeros(6,n);
Xi=transpose(mvnrnd(zero,sigmamatrix,N)); %Initilization

%p = @(x1,x2,y) prod(normpdf(Y,(6*v*ones(1,N)-10*gamma*log10(vecnorm(([x1;x2]-pos_vec(:,1).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,2).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,3).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,4).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,5).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,6).*ones(2,N))))),zeta),2);%1/(2*pi*zeta^2)^3*exp(-1/2*(sum(y)*ones(1,N)-(6*v*ones(1,N)-10*gamma*log10(vecnorm(([x1;x2]-pos_vec(:,1).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,2).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,3).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,4).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,5).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,6).*ones(2,N)))))));
w = log(p(Xi(1,:),Xi(4,:),Y(:,1)));
L=max(w);
w=exp(w-L);
w=w/sum(w);
tau(1,1) = sum(Xi(1,:).*w)/sum(w);
tau(2,1) = sum(Xi(4,:).*w)/sum(w);

Wi=transpose(mvnrnd([0,0],sigma^2*eye(2),N));


Zi=transpose(zn(Commands(:,1),:));
for k = 1:500 % main loop
    Xi=theta*(Xi)+phiz*(Zi)+phiw*(Wi); %Calculate each new state
    w =log(w.*p(Xi(1,:),Xi(4,:),Y(:,k+1))); % weighting REMOVE/INCLUDE LOG
    L=max(w);
    w=exp(w-L);
    w=w/sum(w);
    Zi=transpose(zn(Commands(:,k+1),:));
    Wi=transpose(mvnrnd([0,0],sigma^2*eye(2),N));
    tau(1,k+1) =  sum(Xi(1,:).*w)/sum(w);
    tau(2,k+1) =  sum(Xi(4,:).*w)/sum(w);
end

plot(tau(1,:),tau(2,:))
title("SIS, N=10000")
hold on
plot(pos_vec(1,1),pos_vec(2,1),'d')
hold on
plot(pos_vec(1,2),pos_vec(2,2),'d')
hold on
plot(pos_vec(1,3),pos_vec(2,3),'d')
hold on
plot(pos_vec(1,4),pos_vec(2,4),'d')
hold on
plot(pos_vec(1,5),pos_vec(2,5),'d')
hold on
plot(pos_vec(1,6),pos_vec(2,6),'d')
%figure
%title("Histogram n=70")
%hist(w)
ESSM=1/(sum((w/sum(w)).^2))
%% Part 4


tau=zeros(2,n);

Xi=transpose(mvnrnd(zero,sigmamatrix,N)); %Initilization

%p = @(x1,x2,y) prod(normpdf(Y,(6*v*ones(1,N)-10*gamma*log10(vecnorm(([x1;x2]-pos_vec(:,1).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,2).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,3).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,4).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,5).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,6).*ones(2,N))))),zeta),2);%1/(2*pi*zeta^2)^3*exp(-1/2*(sum(y)*ones(1,N)-(6*v*ones(1,N)-10*gamma*log10(vecnorm(([x1;x2]-pos_vec(:,1).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,2).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,3).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,4).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,5).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,6).*ones(2,N)))))));
w = p(Xi(1,:),Xi(4,:),Y(:,1));
tau(1,1) = sum(Xi(1,:).*w)/sum(w);
tau(2,1) = sum(Xi(4,:).*w)/sum(w);

Wi=transpose(mvnrnd([0,0],sigma^2*eye(2),N));


Zi=transpose(zn(Commands(:,1),:));
for k = 1:500 % main loop
    ind=randsample(N,N,true,w);
    Xi=Xi(:,ind);
    Xi=theta*(Xi)+phiz*(Zi)+phiw*(Wi); %Calculate each new state
    w =(p(Xi(1,:),Xi(4,:),Y(:,k+1))); % weighting REMOVE/INCLUDE LOG
    
    Zi=transpose(zn(Commands(:,k+1),:));
    Wi=transpose(mvnrnd([0,0],sigma^2*eye(2),N));
    tau(1,k+1) =  sum(Xi(1,:).*w)/sum(w);
    tau(2,k+1) =  sum(Xi(4,:).*w)/sum(w);
end

plot(tau(1,:),tau(2,:))
title("SISR N=10000")
hold on
plot(pos_vec(1,1),pos_vec(2,1),'d')
hold on
plot(pos_vec(1,2),pos_vec(2,2),'d')
hold on
plot(pos_vec(1,3),pos_vec(2,3),'d')
hold on
plot(pos_vec(1,4),pos_vec(2,4),'d')
hold on
plot(pos_vec(1,5),pos_vec(2,5),'d')
hold on
plot(pos_vec(1,6),pos_vec(2,6),'d')

%% Part 5

load("RSSI-measurements-unknown-sigma.mat")
T = 100;

tau=zeros(6,n);

Xi=transpose(mvnrnd(zero,sigmamatrix,N)); %Initilization

%p = @(x1,x2,y) prod(normpdf(Y,(6*v*ones(1,N)-10*gamma*log10(vecnorm(([x1;x2]-pos_vec(:,1).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,2).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,3).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,4).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,5).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,6).*ones(2,N))))),zeta),2);%1/(2*pi*zeta^2)^3*exp(-1/2*(sum(y)*ones(1,N)-(6*v*ones(1,N)-10*gamma*log10(vecnorm(([x1;x2]-pos_vec(:,1).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,2).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,3).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,4).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,5).*ones(2,N)))+vecnorm(([x1;x2]-pos_vec(:,6).*ones(2,N)))))));

% choose sigma that gives maximum values
sigmaRes = zeros(n,1);
sigmaV = zeros(T,1);
wmat = zeros(T,N);
for j = 1:T
    sigmaV(j) = rand(1)*3;
    wmat(j,:) = p5(Xi(1,:),Xi(4,:),Y(:,1),sigmaV(j)); % weighting REMOVE/INCLUDE LOG
end
[wmax,wind] = max(sum(wmat,2));
w = wmat(wind,:);
sigmaRes(k+1) = sigmaV(wind);

tau(1,1) = sum(Xi(1,:).*w)/sum(w);
tau(2,1) = sum(Xi(4,:).*w)/sum(w);

Wi=transpose(mvnrnd([0,0],sigma^2*eye(2),N));


Zi=transpose(zn(Commands(:,1),:));
for k = 1:n-1 % main loop
    ind=randsample(N,N,true,w);
    Xi=Xi(:,ind);
    Xi=theta*(Xi)+phiz*(Zi)+phiw*(Wi); %Calculate each new state
    
    % --------- Sigma part ---------
    for j = 1:T
        sigmaV(j) = rand(1)*3;
        wmat(j,:) = (p5(Xi(1,:),Xi(4,:),Y(:,k+1),sigmaV(j))); % weighting REMOVE/INCLUDE LOG
    end
    temp = sum(w,2);
    [wmax,wind] = max(sum(wmat,2));
    w = wmat(wind,:);
    sigmaRes(k+1) = sigmaV(wind);
    % ---------         ----------
    
    Zi=transpose(zn(Commands(:,k+1),:));
    Wi=transpose(mvnrnd([0,0],sigma^2*eye(2),N));
    tau(1,k+1) =  sum(Xi(1,:).*w)/sum(w);
    tau(2,k+1) =  sum(Xi(4,:).*w)/sum(w);
end

plot(tau(1,:),tau(2,:))
hold on
plot(pos_vec(1,1),pos_vec(2,1),'d')
