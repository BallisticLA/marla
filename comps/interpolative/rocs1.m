function[Out] = rocs1(A, k, s, p, axis)
        if axis == 0
            % Row ID
            % Sketch construction stage - alternative options are available in 
            %'../rangefinders'.
            class_A = class(A);
            [~, n] = size(A);
            % By default, a Gaussian random sketching matrix is used.
            % Alternative choices are present in '../../utils/sketching_operators'.
            Omega = randn(n, k + s, class_A);
            [Y, ~] = qr(A * Omega, 0);    
            disp(size(Y));
            for j = 1 : p
                [Y, ~] = qr(A' * Y, 0);
                [Y, ~] = qr(A * Y, 0);
            end
            [~, ~, I] = qr(Y', 'vector');
            Out = I(1 : k);
        elseif axis == 1
            % Column ID
            % Sketch construction stage - alternative options are available in 
            %'../rangefinders'.
            % Lines below represent pass A' into rf1
            class_A = class(A);
            [m, ~] = size(A);
            % By default, a Gaussian random sketching matrix is used.
            % Alternative choices are present in '../../utils/sketching_operators'.
            Omega = randn(m, k + s, class_A);
            [Y, ~] = qr(A' * Omega, 0);    
            for j = 1 : p
                [Y, ~] = qr(A * Y, 0);
                [Y, ~] = qr(A' * Y, 0);
            end
            [~, ~, J] = la.qr(Y, 'vetcor');
            Out = J(1 : k);
        end
end