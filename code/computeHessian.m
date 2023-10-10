function Dr = computeHessian(omega,psi,t)
global P power_pattern h_ul V K N
Dr11=zeros(K,K);
for k=1:K
    Dr11(k,k)=-(((h_ul{k}')*inv(V{k}'*diag(psi)*V{k}+omega(k)*h_ul{k}*h_ul{k}')*h_ul{k})^2+1/(t*omega(k)^2));
end
Dr11=real(Dr11);
Dr22=zeros(N,N);
for n=1:N
    Dr22(n,n)=1/(t*psi(n)^2);
    for m=1:N
        for k=1:K
            vn=V{k}(n,:);
            vn=vn';
            vm=V{k}(m,:);
            vm=vm';
            
            Xi=inv(V{k}'*diag(psi)*V{k});
            Gamma=inv(V{k}'*diag(psi)*V{k}+omega(k)*h_ul{k}*h_ul{k}');
            Dr22(n,m)=Dr22(n,m)-vn'*(Gamma*vm*vm'*Gamma-Xi*vm*vm'*Xi)*vn;
        end
    end
end
Dr22=real(Dr22);
Dr12=zeros(K,N);
for n=1:K
    for m=1:N
        Gamma=inv(V{n}'*diag(psi)*V{n}+omega(n)*h_ul{n}*h_ul{n}');
        vm=V{n}(m,:);
        vm=vm';
        
        Dr12(n,m)=-h_ul{n}'*Gamma*vm*vm'*Gamma*h_ul{n};
    end
end
Dr12=real(Dr12);
Dr=[Dr11          Dr12             -ones(K,1) ;
    Dr12'         Dr22              P*power_pattern;
    ones(1,K)     zeros(1,N)                  0                ];