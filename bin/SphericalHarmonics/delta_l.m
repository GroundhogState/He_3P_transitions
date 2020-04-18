function [dl] = delta_l(L,M)
%SUMMARY: Computes nonvanishing contributions to the 
% summation formula for the spherical harmonics.

dl = []; 

for p = 0:L
    for q = 0:L
        for s = 0:L
            
           c1 = ( (p + q + s) == L );
           c2 = ( (p - q) == M );
           
           if c1 && c2
                dl = [dl;[p q s]];
           end
                
        end
    end
end

    
end
    