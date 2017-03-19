classdef DivGrad < PDE
  methods
    function obj = DivGrad(data, fesH1, fesHDiv)
      div = Op_data_div_id(data.a, fesHDiv, fesH1);
      grad = Op_data_Grad_Id(@(x)1+0*x(:,1), fesH1, fesHDiv);
      idV =  Op_data_Id_Id(@(x)1+0*x(:,1), 0, fesHDiv, fesHDiv);
      f = Fc_Data_Id(data.f, fesH1, 0);
      %
      lhs =  {{}, {div}; ...
             {grad},{idV}};
      rhs = {{f}; {}};
      obj = obj@PDE(lhs, rhs);
    end    
    function R = solve2(obj)
      M = obj.getStiffnessBlock(2,2);
      N = size(M,1);
      %M = spdiags(sum(M,2), 0 , N, N);
      Minv = inv(M);
      Minv(abs(Minv)<1e-12) = 0;
      G = obj.getStiffnessBlock(2,1);
      D = obj.getStiffnessBlock(1,2);
      F = obj.getLoadBlock(1);
      S = obj.getShift(1);
      fd = obj.getFreeDoFsTrial(1);
      A = -D*Minv*G;
      b = F - A*S;
      R = zeros(size(b));
      R(~fd) = S(~fd);
      R(fd) = S(fd) + A(fd,fd)\b(fd);
    end
  end
end