% Function that converts streamfunction to horizontal velocity
function [u] = stream2U(s,Nx,Ny)
    U_top = 1;
    dy = 1/(Ny-1);
    
    % Central difference to obtain horizontal velocity U = d(psi)/dy
    stag_x = Nx;
    stag_y = Ny-1;
    
    % u_stag refers to the velocity in the staggered grid
    for i=1:stag_x
        for j=1:stag_y
            u_stag(i,j) = (s(i,j+1) - s(i,j)) / dy;
        end
    end
    
    % Interpolating staggered velocity back to reference grid
    nny = stag_y - 1;
    for i=1:Nx
        for j=1:nny
            u(i,j) = (u_stag(i,j+1) + u_stag(i,j))/2;
        end
    end
    
    % Adding Boundary conditions for u in reference grid
    u = [zeros(Nx,1) u U_top*ones(Nx,1)];
end
