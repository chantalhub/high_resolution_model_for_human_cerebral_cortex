function [U, S, V] = csvd(u_flow, v_flow)
% =========================================================================
% This function is used to perform a singular vector decomposition of a complex velocity field.
% Inputs:
%       u_flow: the velocity in the x-direction of pixels
%       v_flow: the velocity in the y-direction of pixels 
% Outputs:
%       This function returns numeric unitary matrices U and V with the columns containing the singular vectors, and a diagonal matrix S containing the singular values.
% =========================================================================
    [M, N, P] = size(u_flow);
    valid_indices = ~isnan(u_flow(:, :, 1));
    valid_rows = find(valid_indices);
    non_nan_values_u = reshape(u_flow, [], P);
    non_nan_values_u = non_nan_values_u(valid_rows, :);
    non_nan_values_v = reshape(v_flow, [], P);
    non_nan_values_v = non_nan_values_v(valid_rows, :);
    svdmat = non_nan_values_u + 1i*non_nan_values_v;
    [U, S, V] = svd(svdmat);
end
