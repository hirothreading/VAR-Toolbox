clear;
beep off;

addpath('functions');



load priors
load results
fsize=500; %number of MCMC iterations
HORZ=40;   % horizon
scale=1;   % shock size=scale standard deviations
reps=50;   %number of simulations to estimate IRF
pos=4;  %shock is inserted in equation number pos


irf=zeros(fsize,HORZ,N);

for jgibbs=1:fsize
    tic
    jgibbs
b1=bsave1(jgibbs,:);
b2=bsave2(jgibbs,:);
s1=squeeze(sigmaS1(jgibbs,:,:));
s2=squeeze(sigmaS2(jgibbs,:,:));
tar=tsave(jgibbs,1);
delay=tsave(jgibbs,2);

[Y,X,Ystar]=preparexx( data,L,delay,tarvar,VarBench );
%history
Ym=transform(Y,VarBench.transform);
e1=Ystar<=tar;
e2=Ystar>tar;
T=rows(Y);

A01=chol(s1);
A02=chol(s2);


T1=sum(e1);
T2=sum(e2);
history1=cell(T1,1);
historym1=cell(T1,1);
history2=cell(T2,1);
historym2=cell(T2,1);

LL=max([L,delay,Lx+delay]);
history=cell(T,1);
historym=cell(T,1);

for j=LL:T
    history{j}=Y(j-LL+1:j,:);
    historym{j}=Ym(j-LL+1:j,:);
end


jj1=LL;
jj2=LL;
for j=LL:T
    if e1(j)==1
    history1{jj1}=history{j};
    historym1{jj1}=historym{j};
    jj1=jj1+1;
    else
         history2{jj2}=history{j};
    historym2{jj2}=historym{j};
    jj2=jj2+1;
    end
end


Tx=T2;



ir1t=zeros(Tx-LL+1,HORZ,N);
ir2t=zeros(Tx-LL+1,HORZ,N);
ir3t=zeros(Tx-LL+1,HORZ,N);
ir4t=zeros(Tx-LL+1,HORZ,N);

parfor t=LL:T2
    Y0=history2{t};
    Y0m=historym2{t};
    
    [ ir1,irp1] = getimpulse(N, HORZ, L, Y0, b1, ...
    s1,b2,s2,tar,tarvar,delay,reps,Lx,A01,A02,Y0m,scale,pos,VarBench.transform);
    
ir1t(t,:,:)=(ir1);


end

irf(jgibbs,:,:)=squeeze(mean(ir1t(LL+1:end,:,:)));


end

save('irf1sd_r2','irf');
