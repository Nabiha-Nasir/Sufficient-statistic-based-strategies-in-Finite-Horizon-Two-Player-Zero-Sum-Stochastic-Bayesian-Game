function [col_index_sigma] = sigma_col_index_new(A,B,k,k_present,Ha,t,n_is1)
if t==1
    col_index_sigma=k_present; 
else
        col_index_sigma=sum(n_is1(1:t-1));
        for s=1:t-1            
            ks=Ha(3*s-2);
            as=Ha(3*s-1);
            bs=Ha(3*s);
            col_index_sigma=col_index_sigma+(ks-1)*A*B*k*(A*B*k)^(t-1-s)+(as-1)*B*k*(A*B*k)^(t-1-s)+(bs-1)*k*(A*B*k)^(t-1-s);
        end
            kt=Ha(end);
            col_index_sigma=col_index_sigma+(kt-1)+1;
end
