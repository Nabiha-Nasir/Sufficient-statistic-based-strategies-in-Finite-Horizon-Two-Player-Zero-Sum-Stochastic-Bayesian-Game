function [lset,n_lset] = Lset(T,l)
%%% This function is to get only the state information of player 2 at each
%%% stage
    index=0;
    for l1=1:l
        index=index+1;
        lset{1}(index,:)=[l1];
    end
    n_lset(1)=index;
    for t=2:T
        n_row=size(lset{t-1},1);
             index=0;
             for row=1:n_row % for every history information set
               for l_next=1:l %new info about state
                   index=index+1;
                   lset{t}(index,:)=[lset{t-1}(row,:) l_next];
               end
             end
        n_lset(t)=index; %number of row in each stage
    end
end
