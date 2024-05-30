function mumax = muhat(N,M,nu,nubar)
% muhat(j) computes an upper bound on mu(j,k)
% over all k not in I^+_M.
% It does this for all j in I_N 
% Here mu(j,k) = 1/4*sum_{i=1}^4 1/(nubar^{|j-h_i(k)|}*nu^k)
%
% Input:
% N,M: 1x2 array containing sizes of N,M
% it is assumed M >= N componentwise
% nu: Weight corresponding to the l^1_nu norm (1x2 array), nu>=1
% nubar: Weight used in the aliasing error estimate
% We assume that nubar >= nu
%
% Output:
% mumax: largest value of mu(j,k) over all k of the form
% k=(0:M(1),M(2)+1) and k=(M(1)+1,0:M(2))

% We use that mu(j,k) can be rewritten as 
% 1/4*sum_{i=1}^4 1/(nubar^{|h_i(j)-k|}*nu^k)

[j1,j2] = ndgrid(-N(1):N(1),-N(2):N(2));
inu=1./nu;
inubar=1./nubar;
mumax=intval(zeros(N+1));
musum=intval(zeros([N+1,4]));

invfac2=inu(2)^(M(2)+1)*inubar(2).^(M(2)+1-j2);
for k1=0:M(1)
    invnubarnu=inu(1)^k1*inubar(1).^(abs(k1-j1)).*invfac2;
    musum(:,:,1)=invnubarnu(N(1)+1:2*N(1)+1,N(2)+1:2*N(2)+1);
    musum(:,:,2)=invnubarnu(N(1)+1:-1:1,N(2)+1:2*N(2)+1);
    musum(:,:,3)=invnubarnu(N(1)+1:2*N(1)+1,N(2)+1:-1:1);
    musum(:,:,4)=invnubarnu(N(1)+1:-1:1,N(2)+1:-1:1);
    mumax=max(mumax,sum(musum,3));
end    
invfac1=inu(1)^(M(1)+1)*inubar(1).^(M(1)+1-j1);
for k2=0:M(2)
    invnubarnu=inu(2)^k2*inubar(2).^(abs(k2-j2)).*invfac1;
    musum(:,:,1)=invnubarnu(N(1)+1:2*N(1)+1,N(2)+1:2*N(2)+1);
    musum(:,:,2)=invnubarnu(N(1)+1:-1:1,N(2)+1:2*N(2)+1);
    musum(:,:,3)=invnubarnu(N(1)+1:2*N(1)+1,N(2)+1:-1:1);
    musum(:,:,4)=invnubarnu(N(1)+1:-1:1,N(2)+1:-1:1);
    mumax=max(mumax,sum(musum,3));
end

mumax=mumax/4;

end