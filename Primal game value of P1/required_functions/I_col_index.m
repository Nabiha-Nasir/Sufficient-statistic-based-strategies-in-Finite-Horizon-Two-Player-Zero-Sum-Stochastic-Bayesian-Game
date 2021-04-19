function [col_index_I] = I_col_index(n_is1,A,B,k,I,n_is2,t)
%%% Find the column index for specific information set
    col_index_I=sum(n_is2)*B;
    if t==1
        col_index_I=col_index_I+I; 
    else
        col_index_I=col_index_I+sum(n_is1(1:t-1));
        for s=1:t-1
            ks=I(3*s-2);
            as=I(3*s-1);
            bs=I(3*s);
            col_index_I=col_index_I+(ks-1)*A*B*k*(A*B*k)^(t-1-s)+(as-1)*B*k*(A*B*k)^(t-1-s)+(bs-1)*k*(A*B*k)^(t-1-s);
        end
        kt=I(end);
        col_index_I=col_index_I+(kt-1)+1;
    end
end