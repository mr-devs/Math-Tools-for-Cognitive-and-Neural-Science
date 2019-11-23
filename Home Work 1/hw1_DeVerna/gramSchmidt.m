function [orthonormalMAT] = gramSchmidt(n)
% Follow the gram-schmidt process to create an orthonormal matrix of basis
% vectors

% Create first vector
arb_vec = randn(1,10);

lengthVec = @(vec) sqrt(sum(vec*vec')) ;

% Get length
len_arb_vec = lengthVec(arb_vec) ;

% use that to normalize it
basis1 = arb_vec/len_arb_vec ;


% Make second vector
arb_vec2 = randn(1,10);

% project onto the first basis
firstProjection = (arb_vec2*basis1')*basis1 ;

% Find the orthogonal portion by removing the projection from the first
% vector
orthogonalVec = arb_vec2 - firstProjection ;

% Now take its length...
len_arb_vec2 = lengthVec(orthogonalVec) ;

% Normalize the new vector
basis2 = orthogonalVec/len_arb_vec2 ;

% Check it is orthogonal to the last one
vecLenAngle(basis2,basis1) ;
[len1, len2, Angle] = vecLenAngle(basis2,basis1);

fprintf('\n\nBasis 1 vs. Basis 2\n')

fprintf('vecLen1 = %f \n vecLen2 = %f \n Angle = %f', len1, len2,Angle)


% Make third vector
arb_vec3 = randn(1,10);

% project onto the first basis
secondProjection = (arb_vec3*basis2')*basis2 ;
secondProjection2 = (arb_vec3*basis1')*basis1 ;

% Find the orthogonal portion by removing the projection from the first
% vector
orthogonalVec = arb_vec3 - (secondProjection + secondProjection2) ;

% Now take its length...
len_arb_vec3 = lengthVec(orthogonalVec) ;

% Normalize the new vector
basis3 = orthogonalVec/len_arb_vec3 ;

% Check it is orthogonal to the first one
vecLenAngle(basis1,basis3') ;
[len1, len2, Angle] = vecLenAngle(basis1,basis3') ;

fprintf('\n\nBasis 1 vs. Basis 3\n')

fprintf('vecLen1 = %f \n vecLen2 = %f \n Angle = %f', len1, len2,Angle)

% Check it is orthogonal to the first one
[len1, len2, Angle] = vecLenAngle(basis2,basis3');

fprintf('\n\nBasis 2 vs. Basis 3\n')

fprintf('vecLen1 = %f \n vecLen2 = %f \n Angle = %f', len1, len2,Angle)

orthonormalMAT = [basis1 ; basis2 ; basis3] ;

figure 
hold on

plot3([0 orthonormalMAT(1,1)],[0 orthonormalMAT(1,2)],[0 orthonormalMAT(1,3)])
plot3([0 orthonormalMAT(2,1)],[0 orthonormalMAT(2,2)],[0 orthonormalMAT(2,3)])
plot3([0 orthonormalMAT(3,1)],[0 orthonormalMAT(3,2)],[0 orthonormalMAT(3,3)])


xlabel('x-axis')
ylabel('y-axis')
zlabel('z-axis')
view(3)
axis square
grid on
rotate3d on
legend('basis1','basis2','basis3')

end
