clc;
clear;

%% Input Properties
L = 10000; % mm
E = 200000; % MPa
width = 100; % mm, analysis width
d = 1000; % mm, analysis thickness
t = 10; % mm

displacementMiddle = zeros(1, 50);
for r = 1:1:50
    d = L/r;
theta = 0; % orientation
%I = width*d^3/12; % second moment of area
%A = width*d;
%I = (width*d^3/12) - (((width-t)*(d-(2*t))^3)/12);
%A = (2*t*width) + ((d-(2*t))*t);
I = (width*d^3/12) - ((width-(2*t))*(d-(2*t))^3/12);
A = (width*d) - ((width-(2*t))*(d-(2*t)));
B = 0;
G = 80000; % MPa
%ks = 5/6;
%ks = 0.69;
ks = 0.416;
%% Input Analysis Type
% Choose whether analysis is completed using type 1 = Euler Bernoulli or type 2 = Timoshenko
analysisType = 1;

%% Input the number of elements,nodes and degrees of freedom that will be used
numElements = 100; %elementnum;

% Initialise vector to store displacement values at L/2
%curvatureMiddle = zeros(1, 30);

% Loop through the values of interest
%for r = 1:1:30
%    numElements = r;
    
numNodes = numElements + 1;
numDOF = (4*numElements) + 3; 

%% Input the position of the supports and the support type 
% type1 = pin, type2 = roller, type 3 = fixed
numSupports = 2;
suppPosition = [0 10000]; 
suppType = [1 2];

%% Input distributed loads and position

n = 0;
w = 0;
pLoads = [n;w];

%% Input the position of the loading and the loading type
% type 1 = horizontal load, type 2 = vertical load, type 3 = moment

% initialise counter
%i = 0;
%vertDispMiddle = zeros(1, 19);

loads = [-1000];
%for a = -50000:-25000:-500000
%    loads = a;
loadPosition = [5000];
loadType = [2];

%% Input locations of interest

locations = 0:100:L;
    
    %% Create a list of freedoms for the elements depending on analysis type 

    if analysisType == 1
        % Initialise a vector to store the k matrix and the freedom list for each element
        p = eulerMatrix.empty(0, numElements);

        % Loop through all the elements
        for a = 1:numElements

            % Assign numbers to the degrees of freedom of each element
            freedomList = [4*a - 3, 4*a - 2, 4*a - 1, 4*a, 4*a + 1, 4*a + 2, 4*a + 3];

            % Determine the stiffness and loading matrices for each element and store in p 
            p(a) = eulerMatrix(E, I, A, B, L/numElements, theta, freedomList, pLoads);
        end
    else
        % Initialise a vector to store the k matrix and the freedom list for each element
        p = timoshenkoMatrix.empty(0, numElements);

        % Loop through all the elements
        for a = 1:numElements

            % Assign numbers to the degrees of freedom of each element
            freedomList = [4*a - 3, 4*a - 2, 4*a - 1, 4*a, 4*a + 1, 4*a + 2, 4*a + 3];

            % Determine the stiffness and loading matrices for each element and store in p 
            p(a) = timoshenkoMatrix(E, G, I, A, B, L/numElements, theta, freedomList, pLoads, ks);
        end
    end

    %% Create the global stiffness matrix

    % Initialise empty global stiffness matrix
    K = zeros([(4*numElements) + 3, (4*numElements) + 3]);

    % Loop through every element
    for a = 1:numElements

        % Loop through every dof of each element
        for b = 1:7

            for c = 1:7

                % assign the element stiffness matrix values to the member stiffness matrix

                K(p(a).freedomList(b), p(a).freedomList(c)) = K(p(a).freedomList(b), p(a).freedomList(c)) + p(a).stiffnessMatrix(b,c);

            end
        end
    end

    %% Create a vector to store the member loads

    % Initialise empty loading vector
    Q = zeros((4*numElements) + 3, 1);

    % Loop through every element
    for a = 1:numElements

        % Loop through every dof of each element
        for b = 1:7

            % Assign the element loading vector to the member loading vector
            Q(p(a).freedomList(b), 1) = Q(p(a).freedomList(b),1) + p(a).loadingMatrix(b,1);
        end
    end

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

    %% Partitioning the Stiffness matrix 

    % The K11 matrix is made up of the unrestrained rows and columns of K
    % Create an empty matrix to store K11
    K11 = zeros((numDOF - numRestrained), (numDOF - numRestrained));

    % The K12 matrix is made up of the unrestrained rows and restrained columns of K
    % Create an empty matrix to store K12
    K12 = zeros((numDOF - numRestrained), numRestrained);

    % The K21 matrix is made up of the restrained rows and unrestrained columns of K
    % Create an empty matrix to store K22
    K21 = zeros(numRestrained, (numDOF - numRestrained));

    % The K22 matrix is made up of the restrained rows and columns of K
    % Create an empty matrix to store K22
    K22 = zeros(numRestrained, numRestrained);
        % Create a list of unrestrained nodes
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


    % Begin to partition the matrices
    K11 = K([unrestrainedFreedoms], [unrestrainedFreedoms]);
    K12 = K([unrestrainedFreedoms], [restrainedFreedoms]);
    K21 = K([restrainedFreedoms], [unrestrainedFreedoms]);
    K22 = K([restrainedFreedoms], [restrainedFreedoms]);

    %% Obtain Dk, Qk and Qu

    % All of the restrained DOFs have zero displacement
    displacementK = zeros(numRestrained, 1);

    % Initialise the Du vector to be zeros at all unrestrained DOFs
    displacementU = zeros((numDOF - numRestrained), 1);

    % Initialise the Qk vector to be zeros at all unrestrained DOFs
    forcesK = zeros((numDOF - numRestrained), 1);

    % Initialise the Qu vector to be zeros at all restrained DOFs
    forcesU = zeros(numRestrained, 1);

    % Now need to add in any imposed loads to both Qk and Qu

    % First determine the node and DOF of each load
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

    % Secondly, fill Qk and Qu with the loading values

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
                forcesU(countU) = loads(a);

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
            forcesK(round(loadingDOFs(a) - passedRestrained)) = loads(a);
        end

    end

    % Go through the degrees of freedom and add in the values due to the member
    % loads

    for a = 1:numRestrained
        forcesU(a) = forcesU(a) + Q(restrainedFreedoms(a));
    end

    for a = 1:length(unrestrainedFreedoms)
        forcesK(a) = forcesK(a) + Q(unrestrainedFreedoms(a));
    end

    %% Calculating Du and Qu

    % Calculate Du
    displacementU = K11\(forcesK - (K12 * displacementK));

    % Calculate Qu 
    % Initialise a vector to store the support forces]
    suppForces = zeros(numRestrained, 1);

    suppForces = (K21 * displacementU) + (K22 * displacementK) - forcesU;

    % Create vectors to store all the displacements and forces
    displacement = zeros(1, numDOF);
    force = zeros(1, numDOF);

    %% Constructing displacmement and force vectors 
    % Allocate each displacement/force value to the vectors
    for a = 1:numRestrained
        displacement(restrainedFreedoms(a)) = displacementK(a);
        force(restrainedFreedoms(a)) = suppForces(a);

    end

    for a = 1:(numDOF - numRestrained)
        displacement(unrestrainedFreedoms(a)) = displacementU(a);
        force(unrestrainedFreedoms(a)) = forcesK(a);
    end

    %% Post Processing

    % Create a vector to store the displacements at the locations of interest
    horizDispVector = zeros(1, length(locations));
    vertDispVector = zeros(1, length(locations));
    curvature = zeros(1, length(locations));
    rotVector = zeros(1, length(locations));

    % Create list of node locations
    nodeLocations = zeros(1, numNodes);
    for a = 1:numNodes
        nodeLocations(a) = ((a-1)*L)/(numNodes-1);
    end

    % EULER BERNOULLI
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
            horizDispVector(a) = displacement((4*nodeNumber) - 3);
            vertDispVector(a) = displacement((4*nodeNumber) - 2);
            
            % If we are using Euler-Bernoulli:
            if analysisType == 1

                % Get the DOF number
                if nodeNumber ~= (numElements + 1)
                    vL = displacement(((nodeNumber)*4) - 2);

                    vR = displacement(((nodeNumber+1)*4) - 2);

                    rotL = displacement((nodeNumber*4) - 1);

                    rotR = displacement(((nodeNumber+1)*4) - 1);

                    % Get position of interest as a length between 0 and L 
                    lengthZeroToL = 0;

                    % Get length of element
                    elementL = nodeLocations(nodeNumber+1) - nodeLocations(nodeNumber);
                
                else
                    vL = displacement(((nodeNumber-1)*4) - 2);

                    vR = displacement((nodeNumber*4) - 2);

                    rotL = displacement(((nodeNumber-1)*4) - 1);

                    rotR = displacement((nodeNumber*4) - 1);

                    % Get length of element
                    elementL = nodeLocations(nodeNumber) - nodeLocations(nodeNumber-1);
                    
                    % Get position of interest as a length between 0 and L 
                    lengthZeroToL = elementL;
                end

                curvature(a) = ((((12*lengthZeroToL)/elementL^3)-(6/elementL^2))* vL) + ...
                               (((6*lengthZeroToL/elementL^2) - (4/elementL))* rotL) +...
                               (((6/elementL^2) - (12*lengthZeroToL/elementL^3))* vR) +...
                               (((6*lengthZeroToL/elementL^2) - (2/elementL))* rotR);

            % TIMOSHENKO
            elseif analysisType == 2

                if nodeNumber ~= (numElements + 1)
                    % Get the DOF number for uL, uM etc.

                    rotL = displacement((nodeNumber*4) - 1);

                    rotR = displacement(((nodeNumber+1)*4) - 1);

                    % Get position of interest as a length between 0 and L 
                    lengthZeroToL = 0;

                    % Get length of element
                    elementL = nodeLocations(nodeNumber+1) - nodeLocations(nodeNumber);

                else
                
                    % Get the DOF number for uL, uM etc.

                    rotL = displacement(((nodeNumber-1)*4) - 1);

                    rotR = displacement((nodeNumber*4) - 1);

                    % Get length of element
                    elementL = nodeLocations(nodeNumber) - nodeLocations(nodeNumber-1);
                    
                    % Get position of interest as a length between 0 and L 
                    lengthZeroToL = elementL;
                end
                
                curvature(a) = ((-1 * rotL)/elementL) + (rotR/elementL);
            end
                
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

            % If we are using Euler-Bernoulli:
            if analysisType == 1

                % Get the DOF number for uL, uM etc.

                uL = displacement((leftNode*4) - 3);

                uM = displacement(leftNode*4);

                uR = displacement((rightNode*4) - 3);

                vL = displacement((leftNode*4) - 2);

                vR = displacement((rightNode*4) - 2);

                rotL = displacement((leftNode*4) - 1);

                rotR = displacement((rightNode*4) - 1);

                % Get position of interest as a length between 0 and L 
                lengthZeroToL = (locations(a)- nodeLocations(leftNode));

                % Get length of element
                elementL = nodeLocations(rightNode) - nodeLocations(leftNode);

                % Calculate horizontal and vertical displacement and curvature at the position of
                % interest
                horizDispVector(a) = (((1 - (3*lengthZeroToL/elementL)) + ((2*lengthZeroToL^2)/elementL^2))* uL) + ...
                                     (((4*lengthZeroToL/elementL) - ((4*lengthZeroToL^2)/elementL^2))* uM) +...
                                     (((-lengthZeroToL/elementL) + ((2*lengthZeroToL^2)/elementL^2))* uR);

                vertDispVector(a) = ((((1 - (3*lengthZeroToL^2)/elementL^2) + ((2*lengthZeroToL^3)/elementL^3))* vL) + ...
                                    ((lengthZeroToL - ((2*lengthZeroToL^2)/elementL) + ((lengthZeroToL^3)/elementL^2))* rotL) +...
                                    ((((3*lengthZeroToL^2)/elementL^2) - ((2*lengthZeroToL^3)/elementL^3))* vR) +...
                                    (((-1*(lengthZeroToL^2)/elementL) + ((lengthZeroToL^3)/elementL^2))* rotR));

                curvature(a) = ((((12*lengthZeroToL)/elementL^3)-(6/elementL^2))* vL) + ...
                               (((6*lengthZeroToL/elementL^2) - (4/elementL))* rotL) +...
                               (((6/elementL^2) - (12*lengthZeroToL/elementL^3))* vR) +...
                               (((6*lengthZeroToL/elementL^2) - (2/elementL))* rotR);

            % TIMOSHENKO
            elseif analysisType == 2

                % Get the DOF number for uL, uM etc.

                uL = displacement((leftNode*4) - 3);

                uR = displacement((rightNode*4) - 3);

                vL = displacement((leftNode*4) - 2);

                vM = displacement(leftNode*4);

                vR = displacement((rightNode*4) - 2);

                rotL = displacement((leftNode*4) - 1);

                rotR = displacement((rightNode*4) - 1);

                % Get position of interest as a length between 0 and L 
                lengthZeroToL = (locations(a)- nodeLocations(leftNode));

                % Get length of element
                elementL = nodeLocations(rightNode) - nodeLocations(leftNode);

                % Calculate horizontal and vertical displacement and curvature at the position of
                % interest
                horizDispVector(a) = ((1 - (lengthZeroToL/elementL))* uL) + ...
                                     ((lengthZeroToL/elementL)* uR);

                vertDispVector(a) = ((1 - (3*lengthZeroToL/elementL) + ((2*lengthZeroToL^2)/elementL^2))* vL) + ...
                                    (((4*lengthZeroToL/elementL) - ((4*lengthZeroToL^2)/elementL^2))* vM) +...
                                    (((-1*lengthZeroToL/elementL) + ((2*lengthZeroToL^2)/elementL^2))* vR);

                rotVector(a) = ((1 - (lengthZeroToL/elementL))* rotL) + ...
                               ((lengthZeroToL/elementL)* rotR);
                           
                curvature(a) = ((-1 * rotL)/elementL) + (rotR/elementL);
            end
        end
       
    end

%    verticalDisplacement = displacement(2:4:end);
%    i = i + 1;
    
%    vertDispMiddle(i) = vertDispVector(81);
%end
%    middle = curvature(101)
displacementMiddle(r) = vertDispVector(51);
end 