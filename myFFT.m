%my FFT algorithm for a vector using the Cooley-Tukey approch
function X_k=myFFT(x_n)
j=sqrt(-1);
N=length(x_n); %number of samples

%step 1 binary reversal a_n will be the input to the butterfly
inDe=0:N-1;
inBi=de2bi(inDe,log2(N));
outBi=fliplr(inBi);
outDe=bi2de(outBi)+1;
%create a_n starting array
for i=1:N
    a_n(i)=x_n(outDe(i,1));
end

%step 2 enter the butterfly segmentation
for i=1:log2(N)
    for k=1:power(2,i):N
        for l=0:power(2,i-1)-1
            stage(k+l)=a_n(k+l)+a_n(k+l+power(2,i-1))*exp(-2*pi*j*l/power(2,i));
            stage(k+l+power(2,i-1))=a_n(k+l)+a_n(k+l+power(2,i-1))*-exp(-2*pi*j*l/power(2,i));
        end
    end
    a_n=stage;
end

X_k=a_n;



