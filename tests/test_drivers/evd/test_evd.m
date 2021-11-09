
% QB-backed algorithm
tevd1 = TestEVD1(false);
tevd1.test_fr();
tevd1.test_fp_inexact();
tevd1.test_fp_exact();

% Nystrom, for PSD matrices.
tevd2 = TestEVD2();
tevd2.test_fr();