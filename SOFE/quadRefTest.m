%m = RegularMesh([2,2],[0 1;0 1], 0);
%vm = VisualizerMesh.create(m);
%m.adaptiveRefine(1); m.show(); pause
%m.adaptiveRefine(7); m.show(); pause
%for k = 1:10
%  m.adaptiveRefine(4+k*12); m.show(); pause
%end
%
%m = RegularMesh([2,2],[0 1;0 1], 0);
%vm = VisualizerMesh.create(m);
%for k = 1:15
%  m.adaptiveRefine(rand(m.topology.getNumber(2),1)>0.8);
%end
%figure(1),m.show();
%
N = 30;
sdf = @(x) -(0.25 - sum(abs((x-0.5)./[1 1.2]).^2,2).^(1/2));
m = RegularMesh([N N],[0 1;0 1],0);
%
m = m.removeElementsQuad2Tri(m, sdf);
%m = m.removeElementsTri2Tri(m, sdf);
%
m = m.transformTri2Quad(m);
for k = 1:0
  m.adaptiveRefine(rand(m.topology.getNumber(2),1)>0.8);
end
m.show(); axis([0 1 0 1]); axis tight