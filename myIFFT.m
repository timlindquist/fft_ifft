%my IFFT algorithm for a vector using the inverse approch from myFFT
function x_n=myIFFT(X_k)
j=sqrt(-1);
N=length(X_k); %number of samples

a_n=X_k;

%step 1 reverse through the butterfly segmentation
for i=log2(N):-1:1
    for k=1:power(2,i):N
        for l=0:power(2,i-1)-1
            stage(k+l)=(a_n(k+l)+a_n(k+l+power(2,i-1)))*1/2;
            stage(k+l+power(2,i-1))=(a_n(k+l)-a_n(k+l+power(2,i-1)))*1/(2*exp(-2*pi*j*l/power(2,i)));
        end
    end
    a_n=stage;
end

%step 2 binary reversal from a_n to x_n
inDe=0:N-1;
inBi=de2bi(inDe,log2(N));
outBi=fliplr(inBi);
outDe=bi2de(outBi)+1;
%create a_n starting array
for i=1:N
    x_n(i)=a_n(outDe(i,1));
end


