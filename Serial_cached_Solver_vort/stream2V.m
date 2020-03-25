% Function that converts streamfunction to horizontal velocity
function [v] = stream2V(s, Nx, Ny)
    dx = 1/(Nx-1);
    
    % Central difference to obtain horizontal velocity V = -d(psi)/dx
    stag_x = Nx-1;
    stag_y = Ny;
    
    % u_stag refers to the velocity in the staggered grid
    for i=1:stag_x
        for j=1:stag_y
            v_stag(i,j) = -(s(i+1,j) - s(i,j)) / dx;
        end
    end
    
    % Interpolating staggered velocity back to reference grid
    nnx = stag_x - 1;
    for i=1:nnx
        for j=1:Ny
            v(i,j) = (v_stag(i+1,j) + v_stag(i,j))/2;
        end
    end
    
    % Adding Boundary conditions for u in reference grid
    v = [zeros(1,Ny); v; zeros(1, Ny)];
end
