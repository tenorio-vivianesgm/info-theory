function di=DI_plugin_estimate(x,y,memory)
% function to estimate di between two discrete random processes with
% plug-in method. x and y are stationary, ergodic column vectors with the
% same size (time duration), supposedly jointly Markov with memory no
% greater than "memory"


% building matrices A = Y_{i}, B = Y_{-k}^{-1}, C = X_{-k}^0 
       
    
    N=length(x);
    A=zeros(N-memory,1);
    B=zeros(N-memory,memory);
    C=zeros(N-memory,2*memory+1);
    
    for i=1:N-memory
        A(i,:)=y(i+memory)';
        B(i,:)=y(i:i-1+memory)';
        C(i,:)=[x(i:i+memory); y(i:i-1+memory)]';
    end
   
    di=(JointEntropy_byJuli([A B]) - JointEntropy_byJuli(B)-...
        JointEntropy_byJuli([A C]) + JointEntropy_byJuli(C));
    
   