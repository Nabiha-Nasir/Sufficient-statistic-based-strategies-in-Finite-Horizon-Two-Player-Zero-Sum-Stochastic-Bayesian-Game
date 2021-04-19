function [is2,n_is2]=info_J(T,A,B,l)
%%% This function is to get player 2's information set (is2) and number of total possible information combination at each stage 
    index=0;
    %For 1st stage
    for lj=1:l
        index=index+1;
        is2{1}(index,:)=lj;
    end
    n_is2(1)=index; %number of possible information set at 1st stage 
    %For 2nd to last stage
    for t=2:T
        n_row=size(is2{t-1},1);
        index=0;
        for row=1:n_row % for every possible history information there is a new row
           for a=1:A % to visit each action player player 1 (a)
              for b=1:B % to visit each action player player 2 (b)
                for lj=1:l %to visit each state of player 2
                index=index+1;
                is2{t}(index,:)=[is2{t-1}(row,:) a b lj];
                end
               end
           end
        end
        n_is2(t)=index; %number of possible information set at each stage 
    end
end