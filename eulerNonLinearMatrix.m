classdef eulerNonLinearMatrix
    
    properties
        stiffnessMatrix % stiffness matrix
        freedomList % list of freedoms
        loadingMatrix % loading matrix
        tangentStiffnessMatrix % tangent stiffness matrix
        strainVector
    end
    
    methods % a method is a function within a class
        
        % Create a function that will output the stiffness and mass
        % matrices, and the list of freedoms
        function parameters = eulerNonLinearMatrix(wk, xk, L, numLayers, de, layerYoungsModulus, layerYieldStress, layerYieldStrain, layerYValue, layerArea, theta, freedomList, pLoads)
            
            % Initialise a 7x7 zero matrix to store the stiffness matrix
            parameters.stiffnessMatrix = zeros(7,1);
            
            % Initialise a 1x7 zero matrix to store the loading matrix
            parameters.loadingMatrix = zeros(7,1);
            
            % Initialise a 1x7 zero matrix to store the tangent stiffness matrix
            parameters.tangentStiffnessMatrix = zeros(7,7);
            
            % Initialise b matrix
            bXk = zeros(2,7);
            
            strainGauss = cell(numLayers,1);
            
            % Initialise r matrix, axial force and moment vectors
            r = zeros(2,1);
            axialForce = zeros(1,numLayers);
            moment = zeros(1,numLayers);
            
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
            
            % Initialise vector to store element displacement values
            displacementL = zeros(7,1);
            
            % Extract the displacements for this element
            for a = 1:7
                
                displacementL(a) = de(freedomList(a));
            end
            
            % transform these element displacements from global to local
            % coordinates
            
            displacementL = T * displacementL;
            
            % Initialise a pGauss cell to hold the values associated with
            % each xk
            kGauss = cell(length(xk),1);
    
            % Fill in the stiffness matrix values
            for a = 1:length(xk)
                
                % Fill in the B matrix
                bXk(1,1) = - (3/L) + ((4*xk(a))/L^2);
                bXk(1,4) = (4/L) - ((8*xk(a))/L^2);
                bXk(1,5) = -(1/L) + ((4*xk(a))/L^2);
                bXk(2,2) = ((12*xk(a))/L^3) - (6/L^2);
                bXk(2,3) = ((6*xk(a))/L^2) - (4/L);
                bXk(2,6) = (6/L^2) - ((12*xk(a))/L^3);
                bXk(2,7) = ((6*xk(a))/L^2) - (2/L);
                
                % Calculate N and M
                % Calculate internal actions
                for b = 1:numLayers
                    
                    parameters.strainVector(b) = ([1 (-layerYValue(b))] * bXk * displacementL);
                    strainGauss{b} = [1 (-layerYValue(b))];
                    
                    % is the layer still within the linear elastic range?
                    if abs([1 (-layerYValue(b))] * bXk * displacementL) <= layerYieldStrain(b)

                        % if it is linear elastic, stress = youngs modulus * strain
                        axialForce(b) = layerYoungsModulus(b) * ([1 (-layerYValue(b))] * bXk * displacementL) * layerArea;
                        moment(b) = layerYValue(b) * layerYoungsModulus(b) * ([1 (-layerYValue(b))] * bXk * displacementL) * layerArea;

                    elseif ([1 (-layerYValue(b))] * bXk * displacementL) > layerYieldStrain(b)

                        % if not linear elastic, stress = fy
                        axialForce(b) = layerYieldStress(b) * layerArea;
                        moment(b) = layerYValue(b) * layerYieldStress(b) * layerArea;

                    else 

                        axialForce(b) = -1 * layerYieldStress(b) * layerArea;
                        moment(b) = layerYValue(b) * -1 * layerYieldStress(b) * layerArea;

                    end

                end
                
                % Sum the internal actions at each layer to get the total
                % internal actions and store in r
                
                r(1,1) = sum(axialForce);
                r(2,1) = - sum(moment);
        
                % Calculate the term for this xk value
                kGauss{a} = wk(a) * bXk' * r;
                
                % Add to the total value of qE
                parameters.stiffnessMatrix = parameters.stiffnessMatrix + kGauss{a};
                
            end

            % Calculate the element stiffness matrix - Ke = T' * k * T
            parameters.stiffnessMatrix = T * (L/2 * parameters.stiffnessMatrix);

            % Assign the freedom list to the element
            parameters.freedomList = freedomList;
            
            % Initialise Ne vector
            nE = zeros(2,7);
            
            % Initialise a pGauss cell to hold the values associated with
            % each xk
            qGauss = cell(length(xk),1);
            
            % Calculate the element q matrix
            for a = 1:length(xk)
                
                % Fill in the Ne matrix
                nE(1,1) = 1 - ((3*xk(a))/L) + ((2*xk(a)^2)/L^2);
                nE(1,4) = ((4*xk(a))/L) - ((4*xk(a)^2)/L^2);
                nE(1,5) = -(xk(a)/L) + ((2*xk(a)^2)/L^2);
                nE(2,2) = 1 - ((3*xk(a)^2)/L^2) + ((2*xk(a)^3)/L^3);
                nE(2,3) = xk(a) - ((2*xk(a)^2)/L) + (xk(a)^3/L^2);
                nE(2,6) = ((3*xk(a)^2)/L^2) - ((2*xk(a)^3)/L^3);
                nE(2,7) = - (xk(a)^2/L) + (xk(a)^3/L^2);
                
                % Calculate the term for this xk value
                qGauss{a} = wk(a) * nE' * pLoads;
                
                % Add to the total value of qE
                parameters.loadingMatrix = parameters.loadingMatrix + qGauss{a};
            end
            
            % Calculate the global element loading matrix - Qe = T' * q
            parameters.loadingMatrix = T' * (L/2 * parameters.loadingMatrix);
            
            % Initialise matrix for the derivative of r
            derR = zeros(2,7);
            
            % Initialise a pGauss cell to hold the values associated with
            % each xk
            ktGauss = cell(length(xk),1);
            
            % Initialise vectors to store the layer values of the
            % derivatives
            nd1 = zeros(1,numLayers);
            nd2 = zeros(1,numLayers);
            nd3 = zeros(1,numLayers);
            nd4 = zeros(1,numLayers);
            nd5 = zeros(1,numLayers);
            nd6 = zeros(1,numLayers);
            nd7 = zeros(1,numLayers);
            md1 = zeros(1,numLayers);
            md2 = zeros(1,numLayers);
            md3 = zeros(1,numLayers);
            md4 = zeros(1,numLayers);
            md5 = zeros(1,numLayers);
            md6 = zeros(1,numLayers);
            md7 = zeros(1,numLayers);
            
            % Calculate kt
            % Fill in the stiffness matrix values
            for a = 1:length(xk)
                
                % Fill in the B matrix
                bXk(1,1) = - (3/L) + ((4*xk(a))/L^2);
                bXk(1,4) = (4/L) - ((8*xk(a))/L^2);
                bXk(1,5) = -(1/L) + ((4*xk(a))/L^2);
                bXk(2,2) = ((12*xk(a))/L^3) - (6/L^2);
                bXk(2,3) = ((6*xk(a))/L^2) - (4/L);
                bXk(2,6) = (6/L^2) - ((12*xk(a))/L^3);
                bXk(2,7) = ((6*xk(a))/L^2) - (2/L);
                
                for b = 1:numLayers
        
                    % is the layer still within the linear elastic range?
                    if abs([1 (-layerYValue(b))] * bXk * displacementL) <= layerYieldStrain(b)

                        % if it is linear elastic, stress derivative = E
                        nd1(b) = layerYoungsModulus(b) * layerArea * ((-3/L) + ((4*xk(a))/L^2)); 
                        nd2(b) = layerYoungsModulus(b) * layerArea * (-1 * layerYValue(b)) * (((12*xk(a))/L^3) - (6/L^2));
                        nd3(b) = layerYoungsModulus(b) * layerArea * (-1 * layerYValue(b)) * (((6*xk(a))/L^2) - (4/L));
                        nd4(b) = layerYoungsModulus(b) * layerArea * ((4/L) - ((8*xk(a))/L^2));
                        nd5(b) = layerYoungsModulus(b) * layerArea * ((-1/L) + ((4*xk(a))/L^2));
                        nd6(b) = layerYoungsModulus(b) * layerArea * (-1 * layerYValue(b)) * ((6/L^2) - ((12*xk(a))/L^3));
                        nd7(b) = layerYoungsModulus(b) * layerArea * (-1 * layerYValue(b)) * (((6*xk(a))/L^2) - (2/L));
                        md1(b) = layerYoungsModulus(b) * layerYValue(b) * layerArea * ((-3/L) + ((4*xk(a))/L^2)); 
                        md2(b) = layerYoungsModulus(b) * layerYValue(b) * layerArea * (-1 * layerYValue(b)) * (((12*xk(a))/L^3) - (6/L^2));
                        md3(b) = layerYoungsModulus(b) * layerYValue(b)^2 * layerArea * (((6*xk(a))/L^2) - (4/L));
                        md4(b) = layerYoungsModulus(b) * layerYValue(b) * layerArea * ((4/L) - ((8*xk(a))/L^2));
                        md5(b) = layerYoungsModulus(b) * layerYValue(b) * layerArea * ((-1/L) + ((4*xk(a))/L^2));
                        md6(b) = layerYoungsModulus(b) * layerYValue(b)^2 * layerArea * ((6/L^2) - ((12*xk(a))/L^3));
                        md7(b) = layerYoungsModulus(b) * layerYValue(b)^2 * layerArea * (((6*xk(a))/L^2) - (2/L));
                    
                    else

                        % if not linear elastic, stress derivative = 0
                        
                        nd1(b) = 0;
                        nd2(b) = 0;
                        nd3(b) = 0;
                        nd4(b) = 0;
                        nd5(b) = 0;
                        nd6(b) = 0;
                        nd7(b) = 0;
                        md1(b) = 0;
                        md2(b) = 0;
                        md3(b) = 0;
                        md4(b) = 0;
                        md5(b) = 0;
                        md6(b) = 0;
                        md7(b) = 0;
                    end

                end
                
                % collect the values in derR
                derR(1,1) = sum(nd1);
                derR(1,2) = sum(nd2);
                derR(1,3) = sum(nd3);
                derR(1,4) = sum(nd4);
                derR(1,5) = sum(nd5);
                derR(1,6) = sum(nd6);
                derR(1,7) = sum(nd7);
                derR(2,1) = -1 * sum(md1);
                derR(2,2) = -1 * sum(md2);
                derR(2,3) = sum(md3);
                derR(2,4) = -1 * sum(md4);
                derR(2,5) = -1 * sum(md5);
                derR(2,6) = sum(md6);
                derR(2,7) = sum(md7);
                
                % Calculate the term for this xk value
                ktGauss{a} = wk(a) * bXk' * derR;
                
                % Add to the total value of qE
                parameters.tangentStiffnessMatrix = parameters.tangentStiffnessMatrix + ktGauss{a};
            end
            
            % Calculate the element stiffness matrix - Ket = T' * ket * T
            parameters.tangentStiffnessMatrix = T' * L/2 * parameters.tangentStiffnessMatrix * T;
        end
    end
end
