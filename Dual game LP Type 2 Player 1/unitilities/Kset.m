function [kset,n_kset] = Kset(T,k)
%%% This function is to get only the state information of player 1 at each
%%% stage
    index=0;
    for k1=1:k
        index=index+1;
        kset{1}(index,:)=[k1];
    end
    n_kset(1)=index;
    for t=2:T
        n_row=size(kset{t-1},1);
             index=0;
             for row=1:n_row % for every history information set
               for ki=1:k %new info about state
                   index=index+1;
                   kset{t}(index,:)=[kset{t-1}(row,:) ki];
               end
             end
        n_kset(t)=index; %number of row in each stage
    end
end


