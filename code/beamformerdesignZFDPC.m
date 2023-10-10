clear 
clc
global P power_pattern V  K N h_ul
rng(1)
N=4;
K=4;
P = 10;
P = 10.^(P./10);
power_pattern=ones(N,1);
power_pattern=(power_pattern)/sum(power_pattern); 

H = 1/2*(randn(K,N)+1i*randn(K,N));

V=cell(K,1);
h_eff=cell(K,1);
V{1}=eye(N);
h_eff{1}=H(1,:);
for k=2:K
    V{k}=null(H(1:k-1,:)); 
    h_eff{k}=H(k,:)*V{k};
end
%% 
h_ul=cell(K,1);
for k=1:K
    h_ul{k}=(H(k,:)*V{k})';
end



%% Proposed iteration-point beamformer design
BETA = .5;
ALPHA = .01;
MAXITERS = 100;
Sumrate = zeros(MAXITERS,1);
omega=ones(K,1);psi=ones(N,1);mu=0;
t=1;
mygamma = 50;
n=0;
tic
while(t<1e20)
    for iters =1:MAXITERS
        r = computeresidualerror(omega,psi,mu,t);
        if (norm(r) < 1e-8)
            break;
        end
        Dr=computeHessian(omega,psi,t);
        step = -Dr\r;
        domega = step(1:K);
        dpsi = step(K+(1:N));
        dmu=step(K+N+1:end);
        s = 1;
        newomega = omega+s*domega; newpsi = psi+s*dpsi;
        while ((min(newomega) < 0) || (min(newpsi) < 0)),
            s = BETA*s;
            newomega = omega+s*domega;
            newpsi = psi+s*dpsi;
        end;
        
        newr = computeresidualerror(newomega,newpsi,mu+s*dmu,t);
        while (norm(newr) > (1-ALPHA*s)*norm(r)) % Backtracking line search
            s = BETA*s;
            newomega = omega+s*domega;
            newpsi = psi+s*dpsi;
            newr = computeresidualerror(newomega,newpsi,mu+s*dmu,t);
        end
        omega = omega+s*domega; psi = psi+s*dpsi;
        mu=mu+dmu*s;
        n=n+1;
        Sumrate(n)=0;
        for k=1:K
            Sumrate(n)=Sumrate(n)+log(real(det(V{k}' * diag(psi)*V{k}+h_ul{k}*omega(k)*h_ul{k}')))-log(real(det(V{k}' * diag(psi)*V{k})));
        end

    end
    t=t*mygamma;
end
toc
Sumrate(n:end)=[];
gap=abs(Sumrate-Sumrate(end));
semilogy(1:n-1,gap)
ylabel('Residual error')
xlabel('Iteration')
saveas(gcf,'../results/covergence.png')

%% compare to CVX 
cvx_expert true
tic
cvx_begin  sdp quiet
    variable covarmat(sum([N:-1:N-K+1]),sum([N:-1:N-K+1])) complex hermitian
    obj=0;
    power=0;
    idx=1;
    for k=1:K
        covarmat(idx:idx+N-k,idx:idx+N-k) >= 0;
        obj=obj+log(1+real(h_eff{k} * covarmat(idx:idx+N-k,idx:idx+N-k) * h_eff{k}'));
        power=power+diag(V{k} * covarmat(idx:idx+N-k,idx:idx+N-k) * V{k}');
        idx=idx+N-k+1;
    end
    maximize( (obj));
    subject to
    real(power) <= P*power_pattern;

cvx_end
toc
Sumrate_CVX=real(obj);
abs(Sumrate(end)-Sumrate_CVX)
