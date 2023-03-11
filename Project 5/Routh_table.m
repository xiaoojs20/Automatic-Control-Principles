%% Routh table and stability
%
% This function prints the routh table 
% and analyzes the stability of the system.
%
% Xiao Jinsong
% 2022/10/24

%% Function
function Routhtable=Routh_table(inputVector)
    format long
    %% Complete the first two row
    RouthTableRow = length(inputVector);
    RouthTableColumn = round(RouthTableRow/2);
    RouthTable = zeros(RouthTableRow,RouthTableColumn);
    % Fill the first row of the table
    RouthTable(1,:) = inputVector(1,1:2:RouthTableRow);
    % Fill the second row of the table
    if (rem(RouthTableRow,2) ~= 0)
        RouthTable(2,1:RouthTableColumn - 1) = inputVector(1,2:2:RouthTableRow);
    else
        RouthTable(2,:) = inputVector(1,2:2:RouthTableRow);
    end
    %% calculate Routh table's other rows
    % Set epss as a small value
    miniscule = 0.0001;
    % Calculate other rows of the table
    for i = 3:RouthTableRow   
        % special case: row of all zeros
        if RouthTable(i-1,:) == 0
            order = (RouthTableRow - i);
            orderDEC = 0;
            for j = 1:RouthTableColumn - 1
                RouthTable(i-1,j) = (order - orderDEC) * RouthTable(i-2,j);
                orderDEC = orderDEC + 2;
            end
        end
        
        for j = 1:RouthTableColumn - 1
            % compute each element of the table
            RouthTable(i,j) = ((RouthTable(i-1,1) * RouthTable(i-2,j+1)) - ....
                (RouthTable(i-2,1) * RouthTable(i-1,j+1))) / RouthTable(i-1,1);
        end
        %  special case: zero in the first column
        if RouthTable(i,1) == 0
            RouthTable(i,1) = miniscule;
        end
    end
    
    %% Compute number of right hand side poles(unstable poles)
    unstablePoles = 0;
    % Check change in signs
    for i = 1:RouthTableRow - 1
        if sign(RouthTable(i,1)) * sign(RouthTable(i+1,1)) == -1
            unstablePoles = unstablePoles + 1;
        end
    end

    % Print calculated data on screen
    fprintf("~~~~~~~> Routh Table: <~~~~~~~~~\n");
    Routhtable = RouthTable
    
    % Print the stability result on screen
    if unstablePoles == 0
        fprintf("~~~~~~~> it is a stable system. <~~~~~~~\n");
    else
        fprintf("~~~~~~~> it is an unstable system. <~~~~~~~\n");
    end
    
    fprintf("~~~~~~~> Number of right hand side poles =%2.0f <~~~~~~~\n",unstablePoles);
    sysRoots = roots(inputVector);
    fprintf("~~~~~~~> Roots of the system characteristic equation <~~~~~~~");
    sysRoots
end