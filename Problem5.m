A = diag(2*ones(1,10)) + diag((-1)*ones(1,9),1) + diag((-1)*ones(1,9),-1);
I = eye(10);

v = rand(10, 1);
v = v / norm(v);
shiftv = v;

lambda = 1;
shiftLambda = 1;

tol = 10e-16;

trueEigs = eig(A);
maxEig = trueEigs(end);
secondEig = trueEigs(end-1);
shiftMaxEig = trueEigs(4);
shiftSecondEig = trueEigs(3);

convFactor = abs(secondEig / maxEig);
counter = 0;
shiftConvFactor = abs((1 - shiftMaxEig) / (1 - shiftSecondEig));
shiftCounter = 0;

errors = [];
shiftErrors = [];
convFactors = [];
shiftConvFactors = [];

while (abs(maxEig - lambda) > tol)
    v = A * v;
    v = v / norm(v);

    lambda = v.' * A * v;

    errors(end+1) = abs(maxEig - lambda);

    counter = counter + 1;
    convFactors(end+1) = convFactor^(2 * counter);
end

while (abs(shiftMaxEig - shiftLambda) > tol)
    shiftv = (A - I) \ shiftv;
    shiftv = shiftv / norm(shiftv);

    shiftLambda = shiftv.' * A * shiftv;

    shiftErrors(end+1) = abs(shiftMaxEig - shiftLambda);

    shiftCounter = shiftCounter + 1;
    shiftConvFactors(end+1) = shiftConvFactor^(2 * shiftCounter);
end

eigenVec = v
eigenVal = lambda

shiftEigenVec = shiftv
shiftEigenVal = shiftLambda

numErrors = length(errors);
xvals = 1:numErrors;

shiftNumErrors = length(shiftErrors);
shiftxvals = 1:shiftNumErrors;

%%
semilogy(xvals, errors, "LineWidth", 2)
hold on
semilogy(xvals, convFactors, "Linewidth", 2)
hold off
title("Convergence of the Power Method")
legend("Actual Errors at Each Iteration", "Theoretical Errors at Each Iteration")
xlabel("Number of Iterations")
ylabel("Error")

%%
semilogy(shiftxvals, shiftErrors, "LineWidth", 2)
hold on
semilogy(shiftxvals, shiftConvFactors, "Linewidth", 2)
hold off
title("Convergence of the Inverse Iteration Method")
legend("Actual Errors at Each Iteration", "Theoretical Errors at Each Iteration")
xlabel("Number of Iterations")
ylabel("Error")