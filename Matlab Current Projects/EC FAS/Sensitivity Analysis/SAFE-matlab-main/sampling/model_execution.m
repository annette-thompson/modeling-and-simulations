function [Y1,Y2,Y3,Y4] = model_execution(fun_test,X)

N = size(X,1);
Y1 = nan(N,1) ;
Y2 = nan(N,1) ;
Y3 = nan(N,1) ;
Y4 = nan(N,1) ;

parfor j=1:N %runs in paralel
    fprintf('%i\n',j)
    [Y1(j),Y2(j),Y3(j),Y4(j)] = feval(fun_test,X(j,:));
end

end


