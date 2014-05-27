% Random3Plus1 generate a random configuration for the 3+1 problem
%   
%   Construct a random set of point correspondences and a directional
%   correspondence for the 3+1 problem. 
%   
%   draw - 1 = visualize the scene, 0 = do not visualize
%   nr_points - number of image correspondences to create (3 for 3+1)
%   stddev_pix - pixel error in normalized camera coordinates
%   fov_degrees - camera field of view
%   translation_direction - 0 = random, 1 = [1,0,0]' (to the right), 2 =
%                           [0,0,1]' (forward)
%   min_depth - minimum depth of 3D points
%   max_depth - maximum depth of 3D points
%   stdev_g - error in directional correspondence in normalized camera
%               coordinates
%
%   Uses global scale of the right pose (should be set to 1)
%       global scale 
%
%   Output:
%   Q - image points in the first camera
%   Qp - corresponding image points in the second camera
%   G - direction vector in the first camera
%   Gp - corresponding direction vector in the second camera
%   E - essential matrix 
%   P3D - 3D points from which the correspondences were generated
%   P - first camera matrix (always [I,0]
%   Pp - second camera matrix
%
%   When the scene is plotted, the color scheme is as follows:
%       Green line segments: directional correspondence
%       Red XYZ coordinate frame: first camera
%       Blue XYZ coordinate frame: second camera
%       Yellow segments: image point vectors in first camera
%       Magenta segments: image point vectors in second camera
%       Blue stars: 3D points
%
%   To test the algorithm use:
%       [Q,Qp,G,Gp,E,P3D,P,Pp] = Random3Plus1(1,3);
%       Pstack = Solve3Plus1Action(Q,Qp,G,Gp);
%   The columns of Pstack columns correspond to solutions for Pp in
%   row order.
%
% Author: Oleg Naroditsky
function [Q,Qp,G,Gp,E,P3D,P,Pp] = Random3Plus1(draw,nr_points,stdev_pix,fov_degrees,translation_direction,min_depth,max_depth,stdev_g)

if ~exist('stdev_pix','var')
    % noise parameters:
    stdev_pix = 0.0;
end

if ~exist('translation_direction','var')
    % translation_direction: 0 - random
    translation_direction = 0;
end

if ~exist('fov_degrees','var')
    % camera field of view in degrees
    fov_degrees = 30;
end

if ~exist('min_depth','var')
    % minimum depth of 3D points from first camers
    min_depth = 3;
end

if ~exist('max_depth','var')
    % maximum depth of 3D points from first camers
    max_depth = 10;
end

if ~exist('stdev_g','var')
    % maximum depth of 3D points from first camers
    stdev_g = 0;
end

im_w = 640;

nx = tan(fov_degrees/2*pi/180)*2;

stdev = nx/im_w*stdev_pix;

g = rand(3,1);
g = g/norm(g);

% generate a set of 3d points:
P3D = Make3DPoints(fov_degrees,min_depth,max_depth,nr_points);

Pp = MakeRandomPose(translation_direction);

    
% left pose:
R = eye(3);%orth(rand(3));
t = [0,0,0]';
P  = [R, t];

if 0
% right pose:
in_front = 0; % how many points are in front of the camera:
while in_front ~= nr_points
    Pp = MakeRandomPose(translation_direction);
    Qp = Pp*P3D;
    in_front = sum(Qp(3,:) > 0);
end
end;

% rotation in right cam
Rp = Pp(1:3,1:3);
tp = Pp(1:3,4);

% generate a random rotation with deviation stdev_r:
g_noise = expm(randn*stdev_g*pi/180*skew3((rand(3,1)-0.5)*2));

% direction in left coords
G = R*g;

% direction in the right camera:
Gp = Rp*g;

% contaminate direction with noise:
Gp = g_noise*Gp;

% left image points
Q = P*P3D;
Q = Q./repmat(Q(3,:), 3,1);
Q(1:2,:) = Q(1:2,:)+   randn(2,size(Q,2))*stdev;

% right image points
Qp = Pp*P3D;
Qp = Qp./repmat(Qp(3,:),3,1);
Qp(1:2,:) = Qp(1:2,:)+ randn(2,size(Q,2))*stdev;

Pinv = inv([P;[0,0,0,1]]);
Pinv = Pinv(1:3,:);

Ppinv = inv([Pp;[0,0,0,1]]);
Ppinv = Ppinv(1:3,:);

if draw
    ViewPose(Pinv, 'r');ViewPose(Ppinv,'b');
    PlotVector(G(1:3),Pinv, 'g',0);
    PlotVector(Gp(1:3),Ppinv, 'g',0);
    axis equal;
    hold on;
%    plot3(P3D(1,1:3),P3D(2,1:3),P3D(3,1:3),'b*');
    plot3(P3D(1,1:end),P3D(2,1:end),P3D(3,1:end),'b*');    
    
    PlotVector(Q(:,1),Pinv,'y',1);
    PlotVector(Q(:,2),Pinv,'y',1);
    PlotVector(Q(:,3),Pinv,'y',1);

    PlotVector(Qp(:,1),Ppinv,'m',1);
    PlotVector(Qp(:,2),Ppinv,'m',1);
    PlotVector(Qp(:,3),Ppinv,'m',1);

    % calibration matrix for visualization only
    K = [0.5,0,0.5;0,0.5,0.5;0,0,1];
    cameraBox(1,1,K,P,0.3,'g');
    cameraBox(1,1,K,Pp,0.3,'g');
    
    hold off;
end;

% normalized the translation since it's recovered up to scale:
Pp(1:3,4) = Pp(1:3,4)/norm(Pp(1:3,4));
Rp = Pp(1:3,1:3);
tp = Pp(1:3,4);

% essential matrix
E = skew3(R'*tp)*(R'*Rp);

P = [P;0,0,0,1];
Pp = [Pp;0,0,0,1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% translation_direction: 0 - random, 1 - [1,0,0]' (to the right), 2 = [0,0,1]' (forward)
function Pp = MakeRandomPose(translation_direction)
global scale;
scale = 1;
% rotation in right cam
Rp = orth(rand(3));

Rp = ang2orth((rand(1,3)-0.5).*[100,180,100].*pi/180);

% right pose:
if translation_direction == 0
    tp = rand(3,1)*4-2;
    tp = tp/norm(tp)*scale;
else
    if translation_direction == 1
        tp = -Rp*[1,0,0]'*scale;
    else
        if translation_direction == 2
            tp = -Rp*[0,0,1]'*scale;
        end;
    end;
end;



Pp = [Rp,tp];

function [Qr] = ang2orth(theta) 
% ANG2ORTH generate orthogonal matrix from given angles
Qr = eye(3);
n = 1;
for i = 2 : -1 : 1
   for j = 3 : -1 : i + 1
      t = theta(n); 
      s = sin(t); c = cos(t);
      U = eye(3);
      U(i,i) =  c; U(i,j) = s; 
      U(j,i) = -s; U(j,j) = c;
      Qr = U * Qr; n = n + 1;
   end
end

% generate nr_points random 3D points for a camera with fov in degrees, with depth
% limits min_depth and max_depth
function P3D = Make3DPoints(fov,min_depth,depth,nr_points)
P3D = ones(4,nr_points);

c = tan(fov*pi/180/2);

P3D(3,:) = [min_depth+1+rand(1,nr_points)*(depth-min_depth)];
P3D(1:2,:) = [2*(rand(2,nr_points)-0.5)*c].*[P3D(3,:);P3D(3,:)];


function ViewPose(P, c)

unhold = 0;
if ~ishold
    hold on;
    unhold = true;
end

X = [1 0 0 0 0;
       0 0 1 0 0;
       0 0 0 0 1;
       1 1 1 1 1];
   
Xc = P*X;   

plot3(Xc(1,:), Xc(2,:), Xc(3,:), c);
text(Xc(1,1),Xc(2,1),Xc(3,1), 'x');
text(Xc(1,3),Xc(2,3),Xc(3,3), 'y');
text(Xc(1,5),Xc(2,5),Xc(3,5), 'z');

if unhold
    hold off;
end

function wBounds=cameraBox(w,h,K,P, scale,col)

imBounds = [0,        w,        w     0;
            0,        0,        h,   h;
            1,        1,        1,       1;
            1,        1,        1,       1];
K(1:2,1:2) = K(1:2,1:2)/scale;
P(4,1:4) = [0,0,0,1];
K(4,1:3) = [0,0,0];
K(1:4,4) = [0,0,0,1]';

        
wBounds = inv(P)*inv(K)*imBounds;

width = 2;

line([wBounds(1,1), wBounds(1,2)], [wBounds(2,1), wBounds(2,2)],[wBounds(3,1), wBounds(3,2)],'Color',col,'LineWidth',width);
line([wBounds(1,2), wBounds(1,3)], [wBounds(2,2), wBounds(2,3)],[wBounds(3,2), wBounds(3,3)],'Color',col,'LineWidth',width);
line([wBounds(1,3), wBounds(1,4)], [wBounds(2,3), wBounds(2,4)],[wBounds(3,3), wBounds(3,4)],'Color',col,'LineWidth',width);
line([wBounds(1,4), wBounds(1,1)], [wBounds(2,4), wBounds(2,1)],[wBounds(3,4), wBounds(3,1)],'Color',col,'LineWidth',width);

function PlotVector(V,P,c,normalize)

if normalize
    for i = 1:size(V,2)
        if norm(V(:,i)) ~= 0
            V(:,i) = V(:,i)/norm(V(:,i));
        end;
    end;
end;

v = P*[[V;1],[0,0,0,1]'];

unhold = 0;
if ~ishold
    hold on;
    unhold = true;
end
%arrow3d(v(:,2)',v(:,1)', 15, 'cylinder',[0.1,0.1]);
plot3(v(1,:),v(2,:),v(3,:), c);
if unhold
    hold off;
end

