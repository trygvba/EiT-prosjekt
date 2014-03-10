function U = putDirichletBack(Uprev, lowerPlate, upperPlate, uzilow, uziup)

data = [upperPlate uziup; lowerplate uzilow];
nodes = sort([upperplate; lowerPlate],'ascend');
U = Uprev;
for in=nodes
    i = 3*in;
    temp = data(find(data(:,1))==in,2);
    U = [U(1:(i-1)); temp; U(i:end)];
end


end