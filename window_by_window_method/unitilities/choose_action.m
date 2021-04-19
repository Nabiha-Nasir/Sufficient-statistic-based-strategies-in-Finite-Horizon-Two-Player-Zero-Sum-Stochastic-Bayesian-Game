function [a] = choose_action(sigma,A,col_index_sigma)
%Choose player 1's action a_1 accoriding to the strategy sigma
    rd=rand(1);
    a=1;
     while sum(sigma(1:a,col_index_sigma))<rd
         a=a+1;
         if a>A %For some reason we are getting a>A here. To solve this problem we have used a-1 and break
            a=a-1;
            break
         end
     end
end

