function [col_index_J] = J_col_index(n_is1,A,B,l,J,n_is2,t)
%%% Find the column index for specific information set
    col_index_J=sum(n_is1)*A;
    if t==1
        col_index_J=col_index_J+J; 
    else
        col_index_J=col_index_J+sum(n_is2(1:t-1));
        for s=1:t-1            
                ls=J(3*s-2);
                as=J(3*s-1);
                bs=J(3*s);
                col_index_J=col_index_J+(ls-1)*A*B*l*(A*B*l)^(t-1-s)+(as-1)*B*l*(A*B*l)^(t-1-s)+(bs-1)*l*(A*B*l)^(t-1-s);
        end
        lt=J(end);
        col_index_J=col_index_J+(lt-1)+1;
    end
end
