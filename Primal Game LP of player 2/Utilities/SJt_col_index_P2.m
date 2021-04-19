function [col_index_SJt] = SJt_col_index_P2(t,J,A,B,l,b,n_is2)
%%% Find the column index of S_{J_t} for specific information set J
    if t==1
        col_index_SJt=(J-1)*B+b;
    else
        col_index_SJt=sum(n_is2(1:t-1))*B; %how many columns are there before this stage has started
        for s=1:t-1
            ls=J(3*s-2);
            as=J(3*s-1);
            bs=J(3*s);
            col_index_SJt=col_index_SJt+(ls-1)*A*(B^2)*l*(A*B*l)^(t-1-s)+(as-1)*(B^2)*l*(A*B*l)^(t-1-s)+(bs-1)*B*l*(A*B*l)^(t-1-s);
        end
       lt=J(end);
       col_index_SJt=col_index_SJt+(lt-1)*B+b;
    end
end

