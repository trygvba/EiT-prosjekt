function U = putDirichletBack(Uprev, lowerPlate, upperPlate, uzilow, uziup)

data = [upperPlate uziup; lowerPlate uzilow];
nodes = sort([upperPlate; lowerPlate],'ascend');
U = Uprev;
n=length(nodes);
for j=1:n
    in=nodes(j);
    i = 3*in;
    temp = data(find(data(:,1)==in),2);
    U = [U(1:(i-1)); temp; U(i:end)];
end


end