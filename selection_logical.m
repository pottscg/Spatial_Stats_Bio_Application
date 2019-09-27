%% Process Selection
% Converts Selection from ImageJ into a Logical Array

function Selection_Image = selection_logical(filepath)

% Read .csv file from ImageJ interface

Selection = csvread(filepath,1,0);
X = Selection(:,1);
Y = Selection(:,2);
I = Selection(:,3);

% Remove all selection outside image boundary

idx = (X <= 0) | (Y <= 0);
X(idx) = [];
Y(idx) = [];
I(idx) = [];

% Create sparse representation

Selection_Image = sparse(X,Y,I);

% Rotate to true position

Selection_Image = fliplr(Selection_Image);
Selection_Image = rot90(Selection_Image);

% Convert to Logical

Selection_Image = ~logical(Selection_Image);

%figure; imagesc(Selection_Image);

end