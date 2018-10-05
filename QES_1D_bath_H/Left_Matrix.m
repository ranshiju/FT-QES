function ML=Left_Matrix(M,lambda,D)
Lambda=KRON({diag(lambda),eye(D),diag(lambda)});
ML=Lambda*M;
end