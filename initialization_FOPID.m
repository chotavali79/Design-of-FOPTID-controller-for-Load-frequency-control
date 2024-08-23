function Population=initialization_FOPID(pop,Dim,UB,LB)
Boundary_no= size(UB,1); % numnber of boundaries
% If the boundaries of all variables are equal and user enter a signle
% number for both UB and LB
if Boundary_no==1
     Population=rand(pop,Dim).*(UB-LB)+LB;
end

% If each variable has a different LB and UB
if Boundary_no>1
    for i=1:Dim
        UB_i=UB(i);
        LB_i=LB(i);
        Population(:,i)=rand(pop,1).*(UB_i-LB_i)+LB_i;      
    end
end