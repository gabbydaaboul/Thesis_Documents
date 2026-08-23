classdef eulerMatrix
    
    properties
        stiffnessMatrix % stiffness matrix
        freedomList % list of freedoms
        loadingMatrix % loading matrix
        
    end
    
    methods
        
        function parameters = eulerMatrix(E, I, A, B, L, theta, freedomList, pLoads)
            
            % Initialise a 7x7 zero matrix to store the element stiffness matrix
            parameters.stiffnessMatrix = zeros(7);
            
            % Initialise a 7x1 zero matrix to store the element loading matrix
            parameters.loadingMatrix = zeros(7,1);
            
            % Calculate Ra, Rb, Ri
            RA = E * A;
            RB = E * B;
            RI = E * I;
    
            % Fill in the stiffness matrix values
            parameters.stiffnessMatrix(1, 1) = 7 * RA / (3 * L);
            parameters.stiffnessMatrix(1, 2) = -4 * RB / (L^2);
            parameters.stiffnessMatrix(1, 3) = -3 * RB / L;
            parameters.stiffnessMatrix(1, 4) = -8 * RA / (3 * L);
            parameters.stiffnessMatrix(1, 5) = RA / (3 * L);
            parameters.stiffnessMatrix(1, 6) = 4 * RB / (L^2);
            parameters.stiffnessMatrix(1, 7) = -RB / L;
            parameters.stiffnessMatrix(2, 1) = -4 * RB / (L^2);
            parameters.stiffnessMatrix(2, 2) = 12 * RI / (L^3);
            parameters.stiffnessMatrix(2, 3) = 6 * RI / (L^2);
            parameters.stiffnessMatrix(2, 4) = 8 * RB / (L^2);
            parameters.stiffnessMatrix(2, 5) = -4 * RB / (L^2);
            parameters.stiffnessMatrix(2, 6) = -12 * RI / (L^3);
            parameters.stiffnessMatrix(2, 7) = 6 * RI / (L^2);
            parameters.stiffnessMatrix(3, 1) = -3 * RB / L;
            parameters.stiffnessMatrix(3, 2) = 6 * RI / (L^2);
            parameters.stiffnessMatrix(3, 3) = 4 * RI / L;
            parameters.stiffnessMatrix(3, 4) = 4 * RB / L;
            parameters.stiffnessMatrix(3, 5) = -RB / L;
            parameters.stiffnessMatrix(3, 6) = -6 * RI / (L^2);
            parameters.stiffnessMatrix(3, 7) = 2 * RI / L;
            parameters.stiffnessMatrix(4, 1) = -8 * RA / (3 * L);
            parameters.stiffnessMatrix(4, 2) = 8 * RB / (L^2);
            parameters.stiffnessMatrix(4, 3) = 4 * RB / L;
            parameters.stiffnessMatrix(4, 4) = 16 * RA / (3 * L);
            parameters.stiffnessMatrix(4, 5) = -8 * RA / (3 * L);
            parameters.stiffnessMatrix(4, 6) = -8 * RB / (L^2);
            parameters.stiffnessMatrix(4, 7) = 4 * RB / L;
            parameters.stiffnessMatrix(5, 1) = RA / (3 * L);
            parameters.stiffnessMatrix(5, 2) = -4 * RB / (L^2);
            parameters.stiffnessMatrix(5, 3) = -RB / L;
            parameters.stiffnessMatrix(5, 4) = -8 * RA / (3 * L);
            parameters.stiffnessMatrix(5, 5) = 7 * RA / (3 * L);
            parameters.stiffnessMatrix(5, 6) = 4 * RB / (L^2);
            parameters.stiffnessMatrix(5, 7) = -3 * RB / L;
            parameters.stiffnessMatrix(6, 1) = 4 * RB / (L^2);
            parameters.stiffnessMatrix(6, 2) = -12 * RI / (L^3);
            parameters.stiffnessMatrix(6, 3) = -6 * RI / (L^2);
            parameters.stiffnessMatrix(6, 4) = -8 * RB / (L^2);
            parameters.stiffnessMatrix(6, 5) = 4 * RB / (L^2);
            parameters.stiffnessMatrix(6, 6) = 12 * RI / (L^3);
            parameters.stiffnessMatrix(6, 7) = -6 * RI / (L^2);
            parameters.stiffnessMatrix(7, 1) = -RB / L;
            parameters.stiffnessMatrix(7, 2) = 6 * RI / (L^2);
            parameters.stiffnessMatrix(7, 3) = 2 * RI / L;
            parameters.stiffnessMatrix(7, 4) = 4 * RB / L;
            parameters.stiffnessMatrix(7, 5) = -3 * RB / L;
            parameters.stiffnessMatrix(7, 6) = -6 * RI / (L^2);
            parameters.stiffnessMatrix(7, 7) = 4 * RI / L;
            
            % Initialise a 7x7 zero matrix to store the transformation matrix
            T = zeros(7,7);

            % Fill in the transformation matrix 
            T(1,1) = cos(theta);
            T(1,2) = sin(theta);
            T(2,1) = -sin(theta);
            T(2,2) = cos(theta);
            T(3,3) = 1;
            T(4,4) = 1;
            T(5,5) = cos(theta);
            T(5,6) = sin(theta);
            T(6,5) = -sin(theta);
            T(6,6) = cos(theta);
            T(7,7) = 1;

            % Calculate the element stiffness matrix - Ke = T' * k * T
            parameters.stiffnessMatrix = T' * parameters.stiffnessMatrix * T;

            % Assign the freedom list to the element
            parameters.freedomList = freedomList;
            
            % Calculate the loading vector for this element
            parameters.loadingMatrix(1,1) = L/6 * pLoads(1);
            parameters.loadingMatrix(2,1) = L/2 * pLoads(2);
            parameters.loadingMatrix(3,1) = L^2/12 * pLoads(2);
            parameters.loadingMatrix(4,1) = (2*L)/3 * pLoads(1);
            parameters.loadingMatrix(5,1) = L/6 * pLoads(1);
            parameters.loadingMatrix(6,1) = L/2 * pLoads(2);
            parameters.loadingMatrix(7,1) = (-1*L^2)/12 * pLoads(2);
            
            % Calculate the global element loading matrix - Qe = T' * q
            parameters.loadingMatrix = T' * parameters.loadingMatrix;
        end
    end
end