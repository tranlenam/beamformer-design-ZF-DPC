function r = computeresidualerror(omega,psi,mu,t)
global P power_pattern h_ul V K N
alpha=2;
r3=sum(omega)-P;
r1=zeros(K,1);r2=zeros(N,1);
for k=1:K
    r1(k)=(h_ul{k}')*inv(V{k}'*diag(psi)*V{k}+omega(k)*h_ul{k}*h_ul{k}')*h_ul{k}+1/(t*omega(k));
end
r1=real(r1);
for n=1:N
    r2(n)=-1/(t*psi(n));
    for k=1:K
        v=V{k}(n,:);
        v=v';
        Xi=inv(V{k}'*diag(psi)*V{k});
        r2(n)=r2(n)-v'*omega(k)*Xi*h_ul{k}*h_ul{k}'*Xi*v/(1+omega(k)*h_ul{k}'*Xi*h_ul{k});
    end
end
r2=real(r2);
r=[r1-mu;r2+P*power_pattern*mu;r3];

