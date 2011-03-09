%This is a way to do a circular cross correlation of the columns of an
%array.
%
%In this file we take the FFT of the array, the FFT works on the columns
%only, then we make a second array, which is the conjugate of the first, we
%multiply these 2 arrays in an element by element manner, then do the
%inverse FFT.
%After this we circularly shift the columns of the conjugate array, go
%through the process again, then put the results nextto the previous ones. 
%
%The first results will be autocorrelation, the subsequent ones will be
%circular cross correlation
%
%This goes pretty quickly if the output array is preallocated, so the method
%was developed to allow preallocation.
%
%The input variables are:
%NRows => The number of Rows of the input array
%NColumns=> guess
%Signal => the matrix with the time data in columns, for instance if you
%       are looking at the cross correlation of orthogonal codes, each code
%       stream would be a column, of course this is NRows by NColumns
%
%You can change this m-file to a function if you feel like it

CircXCorrData=zeros(NRows,NColumns^2);%Preallocate space
h = waitbar(0,'Please wait...');
SignalFFT=fft(Signal);
ConjSignalFFT=conj(SignalFFT);
for index= 1 : NColumns : ( (NColumns-1)*NColumns+1 )
    waitbar(index/( (NColumns-1)*NColumns+1 ))
    CircXCorrData( :,index : index+(NColumns-1) )=ifft(SignalFFT.*ConjSignalFFT); 
    ConjSignalFFT=circshift(ConjSignalFFT,[0,1]);
end
close(h)

%Comment the plot out if you want
plot(max(abs(CircXCorrData)))
grid on