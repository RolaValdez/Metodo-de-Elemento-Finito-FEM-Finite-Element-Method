function [M]=transformadalineal(theta)
m={@(x) (cosd(x))^2 @(x) sind(x)*cosd(x) @(x) -(cosd(x))^2 @(x) -sind(x)*cosd(x);...
   @(x) sind(x)*cosd(x) @(x) (sind(x))^2 @(x) -sind(x)*cosd(x) @(x) -(sind(x))^2;...
   @(x) -(cosd(x))^2  @(x) -sind(x)*cosd(x) @(x) (cosd(x))^2 @(x) sind(x)*cosd(x);...
   @(x) -sind(x)*cosd(x) @(x) -(sind(x))^2 @(x) sind(x)*cosd(x) @(x) (sind(x))^2};

M=zeros(4,4);
for i=1:4
    for j=1:4
        M(i,j)=feval(m{i,j},theta);
    end
end
