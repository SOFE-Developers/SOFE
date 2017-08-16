classdef GaussTri < QuadRule
  properties
  end
  methods
    function obj = GaussTri(n)
      obj = obj@QuadRule(n);
    end
    function initData(obj)
      try 
        obj.initDataTab()
      catch
%        warning('Optimal quadrule not supported! Duffy quadrature used');
        %
        [p,w] = obj.getGaussPoints(); 
        p = 0.5 + 0.5*p;
        w = w/2;
        wArray = w*w';
        W = wArray(:);
        P = [kron(ones(length(w),1),p) kron(p, ones(length(w),1))];
        % transform to triangle
        obj.weights = W.*(1-P(:,2));
        obj.points(:,1) = P(:,1).*(1-P(:,2));
        obj.points(:,2) = P(:,2);
      end
    end
    function initDataTab(obj)
      switch obj.order
        case 1
          obj.weights = 1/2;
          obj.points = [1 1]/3;
        case 2
          obj.weights = [1; 1; 1]/6;
          obj.points = [1 1; 1 4; 4 1]/6;
        case 3
          a = 0.2811498024409796;
          b = 0.0521835308923537;
          obj.weights = [a;a;a;b;b;b]/2;
          a = 0.1628828503958919;
          b = 0.4779198835675637;
          obj.points = [a a;1-2*a a; a 1-2*a; ...
                        b b;1-2*b b; b 1-2*b];
        case 4
          a = 0.2233815896780115;
          b = 0.1099517436553219;
          obj.weights = [a;a;a;b;b;b]/2;
          a = 0.4459484909159649;
          b = 0.0915762135097707;
          obj.points = [a a;1-2*a a; a 1-2*a; ...
                        b b;1-2*b b; b 1-2*b];
        case 5
          a = 0.1259391805448272;
          b = 0.1323941527885062;
          obj.weights = [a;a;a;b;b;b;9/40]/2;
          a = 0.1012865073234563;
          b = 0.4701420641051151;
          obj.points = [a a;1-2*a a; a 1-2*a; ...
                        b b;1-2*b b; b 1-2*b; 1/3 1/3];
        case 6
          a = 0.0508449063702068;
          b = 0.1167862757263794;
          c = 0.0828510756183736;
          obj.weights = [a;a;a;b;b;b;c;c;c;c;c;c]/2;
          a = 0.0630890144915022;
          b = 0.2492867451709104;
          c1 = 0.0531450498448169;
          c2 = 0.3103524510337844;
          c3 = 1-c1-c2;
          obj.points = [a a;1-2*a a; a 1-2*a; ...
                        b b;1-2*b b; b 1-2*b; ...
                        c1 c2;c2 c1;c1 c3;c3 c1;c2 c3;c3 c2];
        case 7
          a = 0.0135338625156656;
          b = 0.0789512544320110;
          c = 0.1286079278189061;
          d = 0.0561201442833754;
          obj.weights = [a;a;a;b;b;b;c;c;c;d;d;d;d;d;d]/2;
          a = 0.0282639241560763;
          b = 0.4743113232672226;
          c = 0.2411433258498488;
          c1 = 0.7612227480245238;
          c2 = 0.0462708777988089;
          c3 = 1-c1-c2;
          obj.points = [a a;1-2*a a; a 1-2*a; ...
                        b b;1-2*b b; b 1-2*b; ...
                        c c;1-2*c c; c 1-2*c; ...
                        c1 c2;c2 c1;c1 c3;c3 c1;c2 c3;c3 c2];
        case 8
          z = 0.1443156076777872;
          a = 0.1032173705347183;
          b = 0.0324584976231981;
          c = 0.0950916342672846;
          d = 0.0272303141744350;
          obj.weights = [z;a;a;a;b;b;b;c;c;c;d;d;d;d;d;d]/2;
          a = 0.1705693077517602;
          b = 0.0505472283170310;
          c = 0.4592925882927232;
          c1 = 0.2631128296346381;
          c2 = 0.0083947774099576;
          c3 = 1-c1-c2;
          obj.points = [1/3 1/3;a a;1-2*a a;a 1-2*a; ...
                        b b;1-2*b b;b 1-2*b; ...
                        c c;1-2*c c;c 1-2*c; ...
                        c1 c2;c2 c1;c1 c3;c3 c1;c2 c3;c3 c2];
        otherwise % TODO case 9 ... 29
          error('!Catch me!');
      end
    end
  end
end