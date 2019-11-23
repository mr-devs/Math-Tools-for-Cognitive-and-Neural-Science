function [vec1Len, vec2Len, Angle] = vecLenAngle(vec1,vec2)
% Return the length of two vectors and the angle between them
%   This function takes in any two vectors and then returns the length of each
%   vector. Also, should there be an angle between them (i.e. both vectors
%   have a length > 1) then the angle is returned.

%{
Details: 
    - This function can take row or column vectors.
    - Vectors with zero length are returned as such and it is assumed that
    the angle between the vectors is zero - as you cannot have an angle
    between two vectors when there is only one vector :).

%}

% First make sure that the vectors are both the same length...

if length(vec1) ~= length(vec2)
    error('Vectors must have the same length')
end

% Determine whether it is a column or row vector 
% Change to row vector if a column vector is found...

if size(vec1,1) > size(vec1,2)
    vec1 = reshape(vec1,1,length(vec1));
end

% Do the same thing for vec2
if size(vec2,1) > size(vec2,2)
    vec2 = reshape(vec2,1,length(vec2));
end

length_vec1 = sqrt(vec1*vec1') ;
length_vec2 = sqrt(vec2*vec2');
    
dotOfVecs = vec1*vec2';
Angle = acosd(dotOfVecs/(length_vec1 * length_vec2));

% The below handles the occasional imaginary number that spit out for the angle
% between a zero length vector and a regular vector.

if isreal(Angle) == 0
    Angle = imag(Angle);
end

if isnan(Angle)
    Angle = 0 ;
end

results = table();
vec1Len = length_vec1 ;
vec2Len = length_vec2 ;
Angle = Angle ;

end

