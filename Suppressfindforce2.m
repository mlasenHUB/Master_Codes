function [ w,fs, VMin ] = Suppressfindforce2( we, f, V, s )
%Finds the vector of minimum  MULTIPLE suppression force given the
%weight of each component:
%f2=(alpha1*f1, alpha2*f1, ..., alphaN*f1)
%   alpha: vector of weights
a=f*V;
N=length(a);
%{
%find the min non-zero value in the excitation force
f1=f;
f1((f1(:)==0))=nan;
%make the excitation reference force fs1 less than the min 
r=rand(1);
fs1=r*min(f1);
%}

%create a vector of weights for a suppression force where fs1(the reference)
%is in a random position (w(p)=1)
%weights can range from -1 to 1.
w=-1+2*rand(1,N);
w(randi(N))=1;

%can introduce user defined weights
if norm(we)~=0
    w=we;
end

%kth mode to be suppressed
for k=1:N
    VMin(k)=-a(k)/( dot(w, V(:,k))); %original forces to suppress each mode (f1 in the description)
end

fs1=VMin(s);
fs=w.*fs1; %random force to suppress mode s


end

