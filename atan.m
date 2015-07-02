function f = atan( f )
%TAN   Tangent of a CHEBFUN2.

% Empty check: 
if ( isempty( f ) ) 
    return
end

op = @(x,y) atan( feval( f, x, y ) ); % Resample
f = chebfun2( op, f.domain );        % Call constructor.

end
