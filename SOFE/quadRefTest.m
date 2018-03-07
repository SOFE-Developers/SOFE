N = 30;
sdf = @(x) -(0.25 - sum(abs((x-0.5)./[1 1]).^2,2).^(1/2));
m = RegularMesh([N N],[0 1;0 1],0);
for k = 1:10
  m.adaptiveRefine(rand(m.topology.getNumber(2),1)>0.9);
end
m.show(); axis([0 1 0 1]); axis tight