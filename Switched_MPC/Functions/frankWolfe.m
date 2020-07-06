%Uses the frankWolf method to find the maximum of (x'*Q*x)/(x'*P*x) such
%that the maxima has a norm of 1 and it is in the convex hull of the verticies
%in V. 
%See: http://www.princeton.edu/~yc5/ele522_optimization/lectures/grad_descent_constrained.pdf
function maxVal = frankWolfe(Q, P, V)
gradient  = @(x) ((x'*P*x)*(Q'+Q)*x - (x'*Q*x)*(P'+P)*x)/((x'*P*x)*(x'*P*x));

%Get unit vectors in each corner of wedge
cornerVects = V';
cornerVects = cornerVects(:,vecnorm(cornerVects)>10e-10);
cornerVects = cornerVects./vecnorm(cornerVects);

%Initialize the variables
xNew = cornerVects(:,1);
xCur = zeros(size(cornerVects(:,1)));
itr = 0;

%While the dist between updates is large
while norm(xNew - xCur, 2) > 10e-5
    %Update variable, itr count, and step size
    xCur = xNew;
    itr = itr + 1;
    stepSize = 2/(itr+1);
    
    %Find cornerVect which maximizes the linear approx
    %[~,maxIdx] = max(cornerVects' * gradient(xCur));
    [~,maxIdx] = max(cornerVects' * ((xCur'*P*xCur)*(Q'+Q)*xCur - (xCur'*Q*xCur)*(P'+P)*xCur)/((xCur'*P*xCur)*(xCur'*P*xCur)));
    y = cornerVects(:, maxIdx);
    
    %Take a step toward the max cornerVect 
    xNew = (1-stepSize)*xCur + stepSize*y;
    
    %Normalize the vector
    xNew = xNew/norm(xNew,2);
end
maxVal = (xCur'*Q*xCur)/(xCur'*P*xCur);
end