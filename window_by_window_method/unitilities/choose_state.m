function [k_present] = choose_state(p,k)
%rd=rand(1);
%rng(1,'philox');
rd=rand(1);
    k_present=1;
    while sum(p(1:k_present))<rd
        k_present=k_present+1;
    end
    %The following part in not needed in this code as sum(p(1:k))=1 here 
    if k_present>k
        error('error occurs when generating state k');
    end
end

