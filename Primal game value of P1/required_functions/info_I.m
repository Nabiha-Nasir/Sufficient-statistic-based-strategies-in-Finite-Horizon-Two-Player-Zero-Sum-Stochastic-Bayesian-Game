function [is1,n_is1]=info_I(T,A,B,k)
%%% This function is to get the information set of player 1 (is1) and number of
%%% total possible information combinition of each stage (n_is1)
    index=0;
    % For 1st stage
    for ki=1:k
        index=index+1;
        is1{1}(index,:)=ki;
    end
    n_is1(1)=index;
    % For 2nd to last stage
    for t=2:T
        n_row=size(is1{t-1},1);
        index=0;
        for row=1:n_row  % for every possible history information there is a new row
           for a=1:A     % to visit each action player player 1 (a)
              for b=1:B  % to visit each action player player 2 (b)
                 for ki=1:k %to visit each state of player 1
                   index=index+1;
                   is1{t}(index,:)=[is1{t-1}(row,:) a b ki];
                 end
              end
           end
         end
         n_is1(t)=index; %number of row in each stage 
    end
end

