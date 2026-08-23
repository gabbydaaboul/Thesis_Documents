clc;
clear;

%% Input Properties
L = 8000; % mm
E = 200000; % MPa
width = 100; % mm, analysis width
d = 400; % mm, analysis thickness
theta = 0; % orientation
I = width*d^3/12; % second moment of area
A = width*d;

%% Input the number of elements,nodes and degrees of freedom that will be used
numElements = 100; %elementnum;

% Loop through the values of interest
%vertDispMiddle = zeros(1, 150);

%for r = 1:100
%    numElements = r;
    

numNodes = numElements + 1;
numDOF = (4*numElements) + 3; 

%% Input the position of the supports and the support type 
% type1 = pin, type2 = roller, type 3 = fixed
numSupports = 2;
suppPosition = [0 8000]; 
suppType = [1 2];

%% Input distributed loads and position
n = 0;
w = 0;
pLoads = [n;w];

%% Input the position of the loading and the loading type
% type 1 = horizontal load, type 2 = vertical load, type 3 = moment
loads = [-200000];
loadPosition = [4000];
loadType = [2];

% initialise counter
% i = 0;

%for a = -50000:-25000:-475000
%    loads = a;
%loadPosition = [4000];
%loadType = [2];

%% Create a vector to store the nodes of the supports
suppNodes = zeros(1, numSupports);

% Determine which nodes the supports are located at
for a = 1:numSupports
    
    suppNodes(a) = ((suppPosition(a)/L) * (numNodes - 1)) + 1;
    
end

% Initialise a counter to count the number of restrained degrees of freedom
numRestrained = 0;

% Determine how many degrees of freedom are restrained due to the supports
for a = 1:numSupports
    
    if suppType(a) == 1
        
        numRestrained = numRestrained + 2;

    elseif suppType(a) == 2
        
        numRestrained = numRestrained + 1;
        
    elseif suppType(a) == 3
       
        numRestrained = numRestrained + 3;
        
    end
end

% Create a vector to store the restrained degrees of freedom
restrainedFreedoms = zeros(1, numRestrained);

% Initialise counter to determine position in restrainedFreedoms vector
pos = 1;

% Determine which degrees of freedom are restrained due to the supports
for a = 1:numSupports
    
    % Determine pinned support DOFs
    if suppType(a) == 1
        % List the restrained horizontal DOF
        restrainedFreedoms(pos) = (4 * suppNodes(a)) - 3;
        % List the restrained vertical DOF
        restrainedFreedoms(pos + 1) = (4 * suppNodes(a)) - 2;
        % Move the counter forward two places
        pos = pos + 2;
        
    % Determine roller support DOFs
    elseif suppType(a) == 2
        % List the restrained vertical DOF
        restrainedFreedoms(pos) = (4 * suppNodes(a)) - 2;
        % Move the counter forward one place
        pos = pos + 1;
    
    % Determine fixed support DOFs
    elseif suppType(a) == 3
        % List the restrained horizontal DOF
        restrainedFreedoms(pos) = (4 * suppNodes(a)) - 3;
        % List the restrained vertical DOF
        restrainedFreedoms(pos + 1) = (4 * suppNodes(a)) - 2;
        % List the restrained moment DOF
        restrainedFreedoms(pos + 2) = (4 * suppNodes(a)) - 1;
        % Move the counter forward three places
        pos = pos + 3;
        
    end 
end

%% Determine the node and DOF of each load
% Initialise vectors to hold the nodes and DOFs
loadingNodes = zeros(1, length(loads));
loadingDOFs = zeros(1, length(loads));

for a = 1:length(loads)
    
    loadingNodes(a) = ((loadPosition(a)/L) * (numNodes - 1)) + 1;
    
end

for a = 1:length(loads)
    
    if mod(a,2) ~= 0 
        if loadType(a) == 1
            loadingDOFs(a) = (4 * loadingNodes(a)) - 3;
            
        elseif loadType(a) == 2
            loadingDOFs(a) = (4 * loadingNodes(a)) - 2;
            
        else 
            loadingDOFs(a) = (4 * loadingNodes(a)) - 1;
        end
    else    
        loadingDOFs(a) = (4 * loadingNodes(a));
    end
end

%% Create a list of unrestrained nodes
unrestrainedFreedoms = zeros(1, (numDOF - numRestrained));

% Start a counter to determine the position in the list of unrestrained
% freedoms
count = 1; 

% Loop through all the degrees of freedom 
for a = 1:numDOF
    
    % Reset 'is restrained' check
    isRestrained = false;
    
    % Start a loop to see if the degree of freedom is restrained
    for b = 1:numRestrained
        
        % Check if the degree of freedom is restrained
        if a == restrainedFreedoms(b)
            
            % If it is, list it as a restrained DOF
            isRestrained = true;
            
        end
    end
    
    % If the freedom is not restrained, add it to the list
    if isRestrained == false
        
            % Put that freedom into the list of unrestrained freedoms
            unrestrainedFreedoms(count) = a;
            
            % Continue to fill the unrestrained freedom matrix
            count = count + 1;
    end 
end

%% Layer Variables - Non Linear Analysis

% Number of Layers
numLayers = 80;

%for r = 1:100
%    numLayers = r;
% Layer Thickness
layerThickness = d/numLayers;

% Area of Layers
layerArea = layerThickness * width;

% Material Properties (one entry for each material)
youngsModulusMaterials = [200000];
yieldStressMaterials = [300];

% List what material each layer is made from
layerMaterial = ones(1,numLayers);

% Initialise vectors to store properties of each layer
layerYoungsModulus = zeros(1,numLayers);
layerYieldStress = zeros(1,numLayers);
layerYieldStrain = zeros(1, numLayers);
layerYValue = zeros(1,numLayers);

% layer properties and y values
for b = 1:numLayers
    layerYoungsModulus(b) = youngsModulusMaterials(layerMaterial(b));
    layerYieldStress(b) = yieldStressMaterials(layerMaterial(b));
    layerYieldStrain(b) = layerYieldStress(b)/layerYoungsModulus(b);
    layerYValue(b) = d/2 - layerThickness/2 - ((b-1) * layerThickness);
end

%% Member - Non-Linear Analysis

% Initially assume zero displacement
de = zeros(numDOF,1);

% For three gauss points:
xkLine = [-0.774596669 0 0.774596669];
wk = [0.5555556 0.8888888 0.5555556];

% Initialise vector for xk values
xk = zeros(1, length(xkLine));

% Calculate xk values
for a = 1:length(xkLine)
    xk(a) = (L/(numElements*2)) * (xkLine(a) + 1);
end

%% Begin to iterate

% Initialise counter for the number of iterations and the tolerance level
iteration = 1;
normTol = 1;

% Initialise vector to store previous value of Q
qPrevious = zeros((4*numElements) + 3,1);

while normTol > 0.00001
    
        % Initialise a vector to store the k matrix and the freedom list for each element
        nL = eulerNonLinearMatrix.empty(0, numElements);

        % Loop through all the elements
        for a = 1:numElements

            % Assign numbers to the degrees of freedom of each element
            freedomList = [4*a - 3, 4*a - 2, 4*a - 1, 4*a, 4*a + 1, 4*a + 2, 4*a + 3];
            
            % Determine the stiffness and loading matrices for each element and store in p 
            nL(a) = eulerNonLinearMatrix(wk, xk, L/numElements, numLayers, de, layerYoungsModulus, layerYieldStress, layerYieldStrain, layerYValue, layerArea, theta, freedomList, pLoads);
        end

    %% Assemble stiffness and loading vectors

    % Initialise empty global stiffness and loading vector
    KNonLinear = zeros((4*numElements) + 3,1);
    QNonLinear = zeros((4*numElements) + 3,1);
    Kt = zeros((4*numElements) + 3,(4*numElements)+3);
    
    % Loop through every element
    for a = 1:numElements

        % Loop through every dof of each element
        for b = 1:7

            KNonLinear(nL(a).freedomList(b), 1) = KNonLinear(nL(a).freedomList(b),1) + nL(a).stiffnessMatrix(b,1);
            QNonLinear(nL(a).freedomList(b), 1) = QNonLinear(nL(a).freedomList(b),1) + nL(a).loadingMatrix(b,1);
            
            for c = 1:7
                
                % Assign the element loading vector to the member loading vector
             Kt(nL(a).freedomList(b), nL(a).freedomList(c)) = Kt(nL(a).freedomList(b), nL(a).freedomList(c)) + nL(a).tangentStiffnessMatrix(b,c);
            end
        end
    end
    
%% Calculate Qr and delta D
 
% The Kt11 matrix is made up of the unrestrained rows and columns of K
% Create an empty matrix to store K11
Kt11 = zeros((numDOF - numRestrained), (numDOF - numRestrained));

% The Kt12 matrix is made up of the unrestrained rows and restrained columns of K
% Create an empty matrix to store K12
Kt12 = zeros((numDOF - numRestrained), numRestrained);

% The Kt21 matrix is made up of the restrained rows and unrestrained columns of K
% Create an empty matrix to store K22
Kt21 = zeros(numRestrained, (numDOF - numRestrained));

% The Kt22 matrix is made up of the restrained rows and columns of K
% Create an empty matrix to store K22
Kt22 = zeros(numRestrained, numRestrained);

% Begin to partition the matrices
Kt11 = Kt([unrestrainedFreedoms], [unrestrainedFreedoms]);
Kt12 = Kt([unrestrainedFreedoms], [restrainedFreedoms]);
Kt21 = Kt([restrainedFreedoms], [unrestrainedFreedoms]);
Kt22 = Kt([restrainedFreedoms], [restrainedFreedoms]);

% All of the restrained DOFs have zero change in displacement (still zero)
deltaDK = zeros(numRestrained, 1);

% Initialise the Du vector to be zeros at all unrestrained DOFs
deltaDU = zeros((numDOF - numRestrained), 1);

% Initialise the qrk vector to be zeros at all unrestrained DOFs
qrK = zeros((numDOF - numRestrained), 1);

% Initialise the qru vector to be zeros at all restrained DOFs
qrU = zeros(numRestrained, 1);

% Now need to add in any imposed loads to both Qk and Qu
for a = 1:length(loads)
    
    % Initialise check to see if the degree of freedom is restrained
    isRestrained = false;
    
    % Initialise a counter to determine how many restrained freedoms have
    % passed
    passedRestrained = 0;
    
    % Initialise counters to determine the position in Qu
    countU = 1;
    e = restrainedFreedoms(countU);
    
    % Loop through the restrained degrees of freedom
    for b = 1:numRestrained
        
        % Is the current freedom restrained?
        if loadingDOFs(a) == restrainedFreedoms(b)
            
            % If it is, set it to true
            isRestrained = true;
            
            % Determine position in restrainedFreedoms
            while e ~= restrainedFreedoms(b)
                countU = countU + 1; 
                e = restrainedFreedoms(countU);
            end
            
            % Enter the load into Qu
            qrU(countU) = loads(a);
            
            % Increase the value of the counter by 1
            countU = countU + 1;
            
        end
        
    end
    
    % If the load is at an unrestrained DOF
    if isRestrained == false
        
        % Determine whether there has been a restrained value yet 
        for c = 1:length(restrainedFreedoms)
            
            if loadingDOFs(a) > restrainedFreedoms(c)
                
                passedRestrained = passedRestrained + 1;
                
            end
        end

        % Put the load into the correct position in Qk
        qrK(round(loadingDOFs(a) - passedRestrained)) = loads(a);
    end

end

% Go through the degrees of freedom and add in the values due to the member
% loads and subtract K values

for a = 1:numRestrained
    qrU(a) = qrU(a) + QNonLinear(restrainedFreedoms(a)) - KNonLinear(restrainedFreedoms(a));
end
    
for a = 1:length(unrestrainedFreedoms)
    qrK(a) = qrK(a) + QNonLinear(unrestrainedFreedoms(a)) - KNonLinear(unrestrainedFreedoms(a));
end

% Calculate Du
deltaDU = Kt11\(qrK - (Kt12 * deltaDK));

% Calculate Qu 
% Initialise a vector to store the support forces
restrainedNLForces = zeros(numRestrained, 1);
restrainedNLForces = (Kt21 * deltaDU) + (Kt22 * deltaDK) - qrU;

% Create vectors to store all the displacements and forces
deltaD = zeros(numDOF, 1);
qR = zeros(numDOF,1);

%% Constructing delta d and QR

% Allocate each displacement/force value to the vectors
for a = 1:numRestrained
    deltaD(restrainedFreedoms(a)) = deltaDK(a);
    qR(restrainedFreedoms(a)) = restrainedNLForces(a) + qrU(a);
    
end

for a = 1:(numDOF - numRestrained)
    deltaD(unrestrainedFreedoms(a)) = deltaDU(a);
    qR(unrestrainedFreedoms(a)) = qrK(a);
end

    % Calculate de(i+1)
    de = de + deltaD;
    
    % Save current value of Q
    qPrevious = qR + KNonLinear;
    
    %% Recalculate qR
        % Loop through all the elements
        for a = 1:numElements

            % Assign numbers to the degrees of freedom of each element
            freedomList = [4*a - 3, 4*a - 2, 4*a - 1, 4*a, 4*a + 1, 4*a + 2, 4*a + 3];

            % Determine the stiffness and loading matrices for each element and store in p 
            nL(a) = eulerNonLinearMatrix(wk, xk, L/numElements, numLayers, de, layerYoungsModulus, layerYieldStress, layerYieldStrain, layerYValue, layerArea, theta, freedomList, pLoads);
        end

    % Reset matrices for K, Q and Kt
    KNonLinear = zeros((4*numElements) + 3,1);
    QNonLinear = zeros((4*numElements) + 3,1);
    Kt = zeros((4*numElements) + 3,(4*numElements)+3);
    
    % Loop through every element
    for a = 1:numElements

        % Loop through every dof of each element
        for b = 1:7

            KNonLinear(nL(a).freedomList(b), 1) = KNonLinear(nL(a).freedomList(b),1) + nL(a).stiffnessMatrix(b,1);
            QNonLinear(nL(a).freedomList(b), 1) = QNonLinear(nL(a).freedomList(b),1) + nL(a).loadingMatrix(b,1);
            
            for c = 1:7
                
            % Assign the element loading vector to the member loading vector
            Kt(nL(a).freedomList(b), nL(a).freedomList(c)) = Kt(nL(a).freedomList(b), nL(a).freedomList(c)) + nL(a).tangentStiffnessMatrix(b,c);
            
            end
            
        end
    end
    
%% Calculate Qr(1+1)
 
% The Kt11 matrix is made up of the unrestrained rows and columns of K
% Create an empty matrix to store K11
Kt11 = zeros((numDOF - numRestrained), (numDOF - numRestrained));

% The Kt12 matrix is made up of the unrestrained rows and restrained columns of K
% Create an empty matrix to store K12
Kt12 = zeros((numDOF - numRestrained), numRestrained);

% The Kt21 matrix is made up of the restrained rows and unrestrained columns of K
% Create an empty matrix to store K22
Kt21 = zeros(numRestrained, (numDOF - numRestrained));

% The Kt22 matrix is made up of the restrained rows and columns of K
% Create an empty matrix to store K22
Kt22 = zeros(numRestrained, numRestrained);

% Begin to partition the matrices
Kt11 = Kt([unrestrainedFreedoms], [unrestrainedFreedoms]);
Kt12 = Kt([unrestrainedFreedoms], [restrainedFreedoms]);
Kt21 = Kt([restrainedFreedoms], [unrestrainedFreedoms]);
Kt22 = Kt([restrainedFreedoms], [restrainedFreedoms]);

% All of the restrained DOFs have zero change in displacement (still zero)
deltaDK = zeros(numRestrained, 1);

% Initialise the Du vector to be zeros at all unrestrained DOFs
deltaDU = zeros((numDOF - numRestrained), 1);

% Initialise the qrk vector to be zeros at all unrestrained DOFs
qrK = zeros((numDOF - numRestrained), 1);

% Initialise the qru vector to be zeros at all restrained DOFs
qrU = zeros(numRestrained, 1);

% Now need to add in any imposed loads to both Qk and Qu
for a = 1:length(loads)
    
    % Initialise check to see if the degree of freedom is restrained
    isRestrained = false;
    
    % Initialise a counter to determine how many restrained freedoms have
    % passed
    passedRestrained = 0;
    
    % Initialise counters to determine the position in Qu
    countU = 1;
    e = restrainedFreedoms(countU);
    
    % Loop through the restrained degrees of freedom
    for b = 1:numRestrained
        
        % Is the current freedom restrained?
        if loadingDOFs(a) == restrainedFreedoms(b)
            
            % If it is, set it to true
            isRestrained = true;
            
            % Determine position in restrainedFreedoms
            while e ~= restrainedFreedoms(b)
                countU = countU + 1; 
                e = restrainedFreedoms(countU);
            end
            
            % Enter the load into Qu
            qrU(countU) = loads(a);
            
            % Increase the value of the counter by 1
            countU = countU + 1;
            
        end
        
    end
    
    % If the load is at an unrestrained DOF
    if isRestrained == false
        
        % Determine whether there has been a restrained value yet 
        for c = 1:length(restrainedFreedoms)
            
            if loadingDOFs(a) > restrainedFreedoms(c)
                
                passedRestrained = passedRestrained + 1;
                
            end
        end

        % Put the load into the correct position in Qk
        qrK(round(loadingDOFs(a) - passedRestrained)) = loads(a);
    end

end

% Go through the degrees of freedom and add in the values due to the member
% loads and subtract K values

for a = 1:numRestrained
    qrU(a) = qrU(a) + QNonLinear(restrainedFreedoms(a)) - KNonLinear(restrainedFreedoms(a));
end
    
for a = 1:length(unrestrainedFreedoms)
    qrK(a) = qrK(a) + QNonLinear(unrestrainedFreedoms(a)) - KNonLinear(unrestrainedFreedoms(a));
end

% Calculate Du
deltaDU = Kt11\(qrK - (Kt12 * deltaDK));

% Calculate Qu 
% Initialise a vector to store the support forces
restrainedNLForces = zeros(numRestrained, 1);
restrainedNLForces = (Kt21 * deltaDU) + (Kt22 * deltaDK) - qrU;

%% Construct Qr(i+1)
% Allocate each displacement/force value to the vectors
for a = 1:numRestrained
    qR(restrainedFreedoms(a)) = restrainedNLForces(a) + qrU(a);
    
end
for a = 1:(numDOF - numRestrained)
    qR(unrestrainedFreedoms(a)) = qrK(a);
end

    %% Determine whether another iteration is necessary
    
    % Initialise vectors to hold the sqared Q and qR terms
    squareQ = zeros(1, length(qPrevious));
    squareQr = zeros(1, length(qR));
    
    % go through and square each term in qR
    for a = 1:length(qR)
        squareQr(a) = abs(qR(a))^2;
    end
    
    % go through and square each term in Q
    for a = 1:length(QNonLinear)
        squareQ(a) = abs(qPrevious(a))^2;
    end
        
    % take the sum of all the squared terms
    qrSquared = sum(squareQr);
    qSquared = sum(squareQ);
    
    % Calculate what the normtol value is 
    normTol = sqrt(qrSquared)/sqrt(qSquared);
    
    % raise the iteration counter
    iteration = iteration + 1;
    
end

 %% Post Processing

 % Input locations of interest
    locations = 0:50:L;
    
    % Create a vector to store the displacements at the locations of interest
    curvature = zeros(1, length(locations));
    
    % Create list of node locations
    nodeLocations = zeros(1, numNodes);
    for a = 1:numNodes
        nodeLocations(a) = ((a-1)*L)/(numNodes-1);
    end

    % Loop through each location of interest
    for a = 1:length(locations)

        % Get location
        index = locations(a);

        % check if you are at a node:
        % initiate counter
        atNode = false;

        % loop through node locations
        for b = 1:numNodes

            if index == nodeLocations(b)

                % note that you are at a node
                atNode = true;

                % note what node number you are at
                nodeNumber = b;

            end
        end

        % If at a node:
        if atNode == true

                % Get the DOF number
                if nodeNumber ~= (numElements + 1)
                    vL = de(((nodeNumber)*4) - 2);

                    vR = de(((nodeNumber+1)*4) - 2);

                    rotL = de((nodeNumber*4) - 1);

                    rotR = de(((nodeNumber+1)*4) - 1);

                    % Get position of interest as a length between 0 and L 
                    lengthZeroToL = 0;

                    % Get length of element
                    elementL = nodeLocations(nodeNumber+1) - nodeLocations(nodeNumber);
                
                else
                    vL = de(((nodeNumber-1)*4) - 2);

                    vR = de((nodeNumber*4) - 2);

                    rotL = de(((nodeNumber-1)*4) - 1);

                    rotR = de((nodeNumber*4) - 1);

                    % Get length of element
                    elementL = nodeLocations(nodeNumber) - nodeLocations(nodeNumber-1);
                    
                    % Get position of interest as a length between 0 and L 
                    lengthZeroToL = elementL;
                end

                curvature(a) = ((((12*lengthZeroToL)/elementL^3)-(6/elementL^2))* vL) + ...
                               (((6*lengthZeroToL/elementL^2) - (4/elementL))* rotL) +...
                               (((6/elementL^2) - (12*lengthZeroToL/elementL^3))* vR) +...
                               (((6*lengthZeroToL/elementL^2) - (2/elementL))* rotR);

        % If not at a node:
        else

            % initiate counter to pass through nodeLocations
            whichNode = 1;

            % Determine which nodes the point of interest is between
            while locations(a) >= nodeLocations(whichNode)
                whichNode = whichNode + 1;
            end

            % Set left node number
            leftNode = whichNode - 1;
            % Set right node number
            rightNode = whichNode;

                % Get the DOF number for uL, uM etc.

                uL = de((leftNode*4) - 3);

                uM = de(leftNode*4);

                uR = de((rightNode*4) - 3);

                vL = de((leftNode*4) - 2);

                vR = de((rightNode*4) - 2);

                rotL = de((leftNode*4) - 1);

                rotR = de((rightNode*4) - 1);

                % Get position of interest as a length between 0 and L 
                lengthZeroToL = (locations(a)- nodeLocations(leftNode));

                % Get length of element
                elementL = nodeLocations(rightNode) - nodeLocations(leftNode);

                curvature(a) = ((((12*lengthZeroToL)/elementL^3)-(6/elementL^2))* vL) + ...
                               (((6*lengthZeroToL/elementL^2) - (4/elementL))* rotL) +...
                               (((6/elementL^2) - (12*lengthZeroToL/elementL^3))* vR) +...
                               (((6*lengthZeroToL/elementL^2) - (2/elementL))* rotR);

           
            
        end
       
    end

    verticalDisplacement = de(2:4:end);
%    len = length(verticalDisplacement);
%    if rem(len, 2) == 0
%        indexR = len/2;
%    else 
%        indexR = (len+1)/2;
    %i = i + 1;
%    end
%    vertDispMiddle(r) = verticalDisplacement(indexR);
%    curvatureMiddle(r) = curvature(81);
    
%end