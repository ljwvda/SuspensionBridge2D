function Jnonlinear = nonlinearityMatrix_PowerSeries(a,n,M)
% nonlinearityMatrix_PowerSeries calculates the nonlinear part of DF in
% terms of Fourier coefficients using the method based on power series
% expansion
% Also works for higher dimensions
%
% Input:
% a: Fourier coefficients of dimension n
% n: Dimension
% M: Cut-off value for the power series of the nonlinearity
%
% Output:
% Jnonlinear: Nonlinear part of DF of dimension n^2

if M==1 % Handle case M=1
    Nxy=(n(1)+1)*(n(2)+1);
    [~,~,ajext] = nonlinearity_PowerSeries(a,M,1);
    if isintval(a)
        apad = intval(zeros(Nxy));
    else
        apad = zeros(Nxy);
    end
    apad(1:size(ajext,1),1:size(ajext,2)) = ajext;

    ajext = apad;
else
    [~,~,ajext] = nonlinearity_PowerSeries(a,M,1);
end

[k1,k2] = ndgrid(0:n(1),0:n(2));
Ix=M*n(1)+(1:n(1)+1);
Iy=M*n(2)+(1:n(2)+1);

Jnonlinear=altzeros([n+1,n+1],a(1));

% Convolution
mat0=altzeros(n+1,a(1));

for i=1:length(k1(:))
mat=mat0;
   for jx=[-k1(i) k1(i)]
      for jy=[-k2(i) k2(i)]
          mat=mat+ajext(Ix-jx,Iy-jy); 
      end
   end 

   % Remove duplications (k=0 case)
   if k1(i)==0
      mat=mat/2;
   end
   if k2(i)==0
      mat=mat/2;
   end

   Jnonlinear(:,:,k1(i)+1,k2(i)+1)=mat;
end

end