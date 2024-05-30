function Jnonlinear = nonlinearityMatrix(a,Nfft,Nout)
% nonlinearityMatrix calculates the nonlinear part of DF in terms of
% Fourier coefficients
%
% Input:
% a: Fourier coefficients 
% Nfft: 1x2 array holding the number of rows and columns for this N. It
% determines the size of the Fourier transform of the nonlinear part. The
% Fourier transformations are done in nonlinearityFourier.m, which is 
% called in nonlinearity.m. This outputs only the part we actually need
% Nout: 1x2 array determining the size of the output Jnonlinear
%
% Output:
% Jnonlinear: Nonlinear part of DF of dimension (Nout+1)^4

[~,bbarprime] = nonlinearity(a,1,Nfft,2*Nout); % Output is of size 4Nout+1

% Tools for convolution
[k1,k2] = ndgrid(0:Nout(1),0:Nout(2));
Ix=2*Nout(1)+1+(0:Nout(1)); 
Iy=2*Nout(2)+1+(0:Nout(2));

Jnonlinear=altzeros([Nout+1,Nout+1],a(1));

% Convolution
mat0=altzeros(Nout+1,a(1)); % Initialize

for i=1:length(k1(:))
    mat=mat0;
    for jx=[-k1(i) k1(i)]
        for jy=[-k2(i) k2(i)]
            mat=mat+bbarprime(Ix-jx,Iy-jy);% Ix-jx is in [1,3Nout+1]
        end
    end
    
    % Remove duplications (k=0 case)
    if k1(i)==0 && k2(i)==0
        mat=mat/4;
    elseif k1(i)==0 || k2(i)==0
        mat=mat/2;
    end
    
    Jnonlinear(:,:,k1(i)+1,k2(i)+1)=mat;
end

end