function random_rows = jitter_sample(scan_dim, grid_size)
%Generates jittered sampling indices. A two dimensional space defined by
%scan_dim is divided into a grid pattern defined by grid_size. One location
%index is randomly chosen from each grid cell. 
%
%Inputs - scan_dim     : Actual dimensions of the spatial scan (DEMO data:
%                        100 by 100 points) 
%         grid_size    : The number of division in the x and y directions 
%                        [x,y] 
%Outputs - random_rows : Indices indicating the random location chosen in
%                        each grid cell. These indices indicate the rows
%                        that would fall into those locations if the two
%                        dimensional scan were flatened into a one
%                        dimension. As such, we are able to use these
%                        "random_rows" indices to mimic a jittered sampling
%                        for one dimensional sampling.     

%CALCULATE JITTERED SAMPLE INDICES
x_grid=ceil(linspace(1,scan_dim(1),grid_size(1)+1));
y_grid=ceil(linspace(1,scan_dim(2),grid_size(2)+1));

x_points=zeros(grid_size);
y_points=zeros(grid_size);
random_rows=zeros(1,grid_size(1)*grid_size(2));
for i = 1:grid_size(1) %For each x_grid
    for j = 1:grid_size(2) %For each y_grid
        x_rand=rand(1);
        y_rand=rand(1);
        
        x=round((x_rand*(x_grid(i+1)-x_grid(i)))+x_grid(i));
        y=round((y_rand*(y_grid(j+1)-y_grid(j)))+y_grid(j));
        
        x_points(i,j)=x;
        y_points(i,j)=y;
        
        if mod(y,2) == 0
            random_rows(grid_size(2)*(i-1)+j)=scan_dim(1)*(y)-x+1;
        else
            random_rows(grid_size(2)*(i-1)+j)=scan_dim(1)*(y-1)+x;
        end
    end
end