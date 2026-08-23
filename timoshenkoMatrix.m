classdef timoshenkoMatrix
    
    properties
        stiffnessMatrix % stiffness matrix
        freedomList % list of freedoms
        loadingMatrix % loading matrix
        
    end
    
    methods % a method is a function within a class
        
        % Create a function that will output the stiffness and mass
        % matrices, and the list of freedoms
        function parameters = timoshenkoMatrix(E, G, I, A, B, L, theta, freedomList, pLoads, ks)
            
            % Initialise a 7x7 zero matrix to store the stiffness matrix
            parameters.stiffnessMatrix = zeros(7);
            
            % Initialise a 7x1 zero matrix to store the loading matrix
            parameters.loadingMatrix = zeros(7,1);
            
            % Calculate Ra, Rb, Ri
            RA = E * A;
            RB = E * B;
            RI = E * I;
            RS = ks * A * G;
    
            % Fill in the stiffness matrix values
            parameters.stiffnessMatrix(1, 1) = RA / L;
            parameters.stiffnessMatrix(1, 3) = -RB / L;
            parameters.stiffnessMatrix(1, 5) = -RA / L;
            parameters.stiffnessMatrix(1, 7) = RB / L;
            parameters.stiffnessMatrix(2, 2) = 7 * RS / (3 * L);
            parameters.stiffnessMatrix(2, 3) = 5 * RS / 6;
            parameters.stiffnessMatrix(2, 4) = -8 * RS / (3 * L);
            parameters.stiffnessMatrix(2, 6) = RS / (3 * L);
            parameters.stiffnessMatrix(2, 7) = RS / 6;
            parameters.stiffnessMatrix(3, 1) = -RB / L;
            parameters.stiffnessMatrix(3, 2) = 5 * RS / 6;
            parameters.stiffnessMatrix(3, 3) = (L * RS / 3) + (RI / L);
            parameters.stiffnessMatrix(3, 4) = -2 * RS / 3;
            parameters.stiffnessMatrix(3, 5) = RB / L;
            parameters.stiffnessMatrix(3, 6) = -RS / 6;
            parameters.stiffnessMatrix(3, 7) = (L * RS / 6) - (RI / L);
            parameters.stiffnessMatrix(4, 2) = -8 * RS / (3 * L);
            parameters.stiffnessMatrix(4, 3) = -2 * RS / 3;
            parameters.stiffnessMatrix(4, 4) = 16 * RS / (3 * L);
            parameters.stiffnessMatrix(4, 6) = -8 * RS / (3* L);
            parameters.stiffnessMatrix(4, 7) = 2 * RS / 3;
            parameters.stiffnessMatrix(5, 1) = -RA / L;
            parameters.stiffnessMatrix(5, 3) = RB / L;
            parameters.stiffnessMatrix(5, 5) = RA / L;
            parameters.stiffnessMatrix(5, 7) = -RB / L;
            parameters.stiffnessMatrix(6, 2) = RS / (3 * L);
            parameters.stiffnessMatrix(6, 3) = -RS / 6;
            parameters.stiffnessMatrix(6, 4) = -8 * RS / (3 * L);
            parameters.stiffnessMatrix(6, 6) = 7 * RS / (3 * L);
            parameters.stiffnessMatrix(6, 7) = -5 * RS / 6;
            parameters.stiffnessMatrix(7, 1) = RB / L;
            parameters.stiffnessMatrix(7, 2) = RS / 6;
            parameters.stiffnessMatrix(7, 3) = (L * RS / 6) - (RI / L);
            parameters.stiffnessMatrix(7, 4) = 2 * RS / 3;
            parameters.stiffnessMatrix(7, 5) = -RB / L;
            parameters.stiffnessMatrix(7, 6) = -5 * RS / 6;
            parameters.stiffnessMatrix(7, 7) = (L * RS / 3) + (RI / L);
            
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
            parameters.loadingMatrix(1,1) = L/2 * pLoads(1);
            parameters.loadingMatrix(2,1) = L/6 * pLoads(2);
            parameters.loadingMatrix(4,1) = (2*L)/3 * pLoads(2);
            parameters.loadingMatrix(5,1) = L/2 * pLoads(1);
            parameters.loadingMatrix(6,1) = L/6 * pLoads(2);
            
            % Calculate the global element loading matrix - Qe = T' * q
            parameters.loadingMatrix = T' * parameters.loadingMatrix;
        end
    end
end