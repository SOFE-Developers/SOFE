classdef MeshTopologyQuadAdapt < MeshTopologyQuad
  properties
    faceStat
  end
  methods % constructor
    function obj = MeshTopologyQuadAdapt(elem)
      obj = obj@MeshTopologyQuad(elem);
      obj.faceStat = zeros(obj.getNumber(1),2);;
    end
  end
  methods % refinement
    function R = adaptiveRefine(obj, idxE)
      if numel(idxE)<obj.getNumber(2) || max(idxE)>1
        idxE = accumarray(idxE,1, [obj.getNumber(2) 1])>0;
      end
      e2F = obj.connectivity{3,2};
      while 1 % 1) mark faces until closure
        idxF = accumarray(reshape(e2F(idxE,:),[],1),1,[obj.getNumber(1) 1]);
        isCh = idxF&obj.faceStat(:,1)<0;
        father = obj.faceStat(isCh,2);
        if all(idxF(father(father>0))), break; end
        idxF(obj.faceStat(isCh,2)) = 1;
        f2E = obj.getFace2Elem(); % TODO: more efficient/local alternative???
        idxE(sum(f2E(obj.faceStat(isCh,2),:),2)) = 1;
      end
      I = idxF & (obj.isBoundary() | (obj.faceStat(:,1)>0));
      idxF(I) = -idxF(I); % removable faces are negative
      isIrr = idxF&obj.faceStat(:,1)>0;
      inIrr = idxF&obj.faceStat(:,1)<=0;
      % 2) node projector
      el = obj.getEntity('0'); fc = obj.getEntity('1'); nN = obj.getNumber(0);
      nNnewF = sum(inIrr); nNnewE = sum(idxE);
      R = [speye(nN); sparse(repmat((1:nNnewF)',1,2), fc(inIrr,:), 0.5, nNnewF, nN); ...
                      sparse(repmat((1:nNnewE)',1,4), el(idxE,:), 0.25, nNnewE, nN)];
      % 3) assign node indices
      idxN = zeros(size(idxF));
      idxN(isIrr) = obj.faceStat(isIrr,1);
      idxN(inIrr) = nN + (1:nNnewF)';
      idxNE = nN+nNnewF+(1:nNnewE)';
      % 4) REFINE:
      % 4a) elem
      E = [el(idxE,:), reshape(idxN(e2F(idxE,:)),[],4), idxNE];
      el(idxE,:) = E(:,[1 5 7 9]);
      el = [el; E(:,[5 2 9 8]); E(:,[7 9 3 6]); E(:,[9 8 6 4])];
      % 4b) face
      nF0 = size(fc,1);
      tmp(:,:,2) = [fc(inIrr,2) idxN(inIrr)];
      tmp(:,:,1) = [fc(inIrr,1) idxN(inIrr)]; % nFx2xnCh --> (nCh*nF)x2
      fc = [fc; reshape(permute(tmp, [2 3 1]), 2, [])'; ...
               [idxN(e2F(idxE,1)) idxNE]; [idxN(e2F(idxE,2)) idxNE]; ...
               [idxN(e2F(idxE,3)) idxNE]; [idxN(e2F(idxE,4)) idxNE]];
      % 4c) e2F
      nF1 = sum(inIrr); nF2 = sum(idxE); rangeF = nF0+2*nF1 + (1:nF2)';
      oo = obj.getOrientation()>0;
      e2FStat = reshape(obj.faceStat(e2F,2), size(e2F));
      inIrr = cumsum(inIrr); e2InIrr = inIrr(e2F);
      I1 = idxE&isIrr(e2F); I2 = idxE&e2InIrr;
      % elem1
      e2F(I2(:,1),1) = nF0 + 2*e2InIrr(I2(:,1),1) - oo(I2(:,1),1);      
      e2F(I1(:,1),1) = e2FStat(I1(:,1),1) + (~oo(I1(:,1),1));
      e2F(idxE,2) = rangeF+2*nF2;
      e2F(I2(:,3),3) = nF0 + 2*e2InIrr(I2(:,3),3) - oo(I2(:,3),3);
      e2F(I1(:,3),3) = e2FStat(I1(:,3),3) + (~oo(I1(:,3),3));
      e2F(idxE,4) = rangeF;
      % elem2
      tmp = cell(3,1);
      for k = 1:3
        tmp{k} = zeros(sum(idxE),4);
      end
      tmp{1}(I2(idxE,1),1) = nF0 + 2*e2InIrr(I2(:,1),1) - (~oo(I2(:,1),1));
      tmp{1}(I1(idxE,1),1) = e2FStat(I1(:,1),1) + oo(I1(:,1),1);
      tmp{1}(:,2) = rangeF+3*nF2;
      tmp{1}(:,3) = rangeF;
      tmp{1}(I2(idxE,4),4) = nF0 + 2*e2InIrr(I2(:,4),4) - oo(I2(:,4),4);
      tmp{1}(I1(idxE,4),4) = e2FStat(I1(:,4),4) + (~oo(I1(:,4),4));
      % elem3
      tmp{2}(:,1) = rangeF+2*nF2;
      tmp{2}(I2(idxE,2),2) = nF0 + 2*e2InIrr(I2(:,2),2) - oo(I2(:,2),2);
      tmp{2}(I1(idxE,2),2) = e2FStat(I1(:,2),2) + (~oo(I1(:,2),2));
      tmp{2}(I2(idxE,3),3) = nF0 + 2*e2InIrr(I2(:,3),3) - (~oo(I2(:,3),3));
      tmp{2}(I1(idxE,3),3) = e2FStat(I1(:,3),3) + oo(I1(:,3),3);
      tmp{2}(:,4) = rangeF+nF2;
      % elem4
      tmp{3}(:,1) = rangeF+3*nF2;
      tmp{3}(I2(idxE,2),2) = nF0 + 2*e2InIrr(I2(:,2),2) - (~oo(I2(:,2),2));
      tmp{3}(I1(idxE,2),2) = e2FStat(I1(:,2),2) + oo(I1(:,2),2);
      tmp{3}(:,3) = rangeF+nF2;
      tmp{3}(I2(idxE,4),4) = nF0 + 2*e2InIrr(I2(:,4),4) - (~oo(I2(:,4),4));
      tmp{3}(I1(idxE,4),4) = e2FStat(I1(:,4),4) + oo(I1(:,4),4);
      e2F = [e2F; cell2mat(tmp)];
      % 4d) faceStat
      idxF1 = (idxF==1);
      obj.faceStat(idxF1,:) = [idxN(idxF1), nF0+2*inIrr(idxF1)-1];
      rmF = (idxF<0) | (idxF==2);
      rmCh = obj.faceStat(rmF,2); rmCh = rmCh(rmCh>0); rmCh = [rmCh; rmCh+1];
      rmCh = rmCh(idxF(rmCh)==0);
      obj.faceStat(rmCh,:) = 0;
      tmp1 = 2*inIrr(idxF1); tmp2 = find(idxF1);
      fStatNew = zeros(2*nF1,2);
      fStatNew(tmp1-1,1) = -1;
      fStatNew(tmp1-1,2) = tmp2;
      fStatNew(tmp1,  1) = -2;
      fStatNew(tmp1,  2) = tmp2;
      obj.faceStat = [obj.faceStat; fStatNew; zeros(4*nF2,2)];
      % 5) remove faces
      cnt = [cumsum(~rmF);nF0-sum(rmF)+(1:2*nF1+4*nF2)'];
      e2F = cnt(e2F);
      I = obj.faceStat(:,2)>0;
      obj.faceStat(I,2) = cnt(obj.faceStat(I,2));
      fc(rmF,:) = []; obj.faceStat(rmF,:) = [];
      % 6) update
      obj.connectivity{obj.dimP+1,1} = el;
      obj.connectivity{2,1} = fc;
      obj.connectivity{3,2} = e2F;
      obj.connectivity{1,1} = (1:max(fc(:)))';
      obj.connectivity{2,2} = (1:size(fc,1))';
      obj.connectivity{3,3} = (1:size(el,1))';
      obj.connectivity{2,3} = [];
    end
  end
  methods % mesh information
    function R = isBoundary(obj)
      R = isBoundary@MeshTopologyQuad(obj);
      R = R & (obj.faceStat(:,1)==0);
    end
  end
end