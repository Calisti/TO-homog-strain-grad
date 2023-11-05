%**************************************************************************
% Categorises boundary/interior nodes & dofs of microcell mesh
%**************************************************************************
% DESCRIPTION
% Splits the nodes and degrees of freedom of a given rectangular micro-
% cell mesh into sets of nodes and degrees of freedom associated with
% interior, boundary and different ("plus", "minus" and corner) portions
% of the boundary.
% It also checks for one-to-one correspondance between nodes on opposing
% sides of the micro-cell rectangle, which is required in the present
% computational treatment of models under the assumption of periodic
% boundary displacement fluctuations with anti-periodic tractions.
% 
% INPUT
%  - coord:  array of microcell nodal coordinates
%  - ndofn:  number of degrees of freedom per node
%  - npoin:  total number of nodes in microcell mesh
%  - option: option for prescribing displacements
%
% OUTPUT
%  - cell_volume: area of rectangular microcell
%  - dofs:       data structure containing arrays of dofs in each category
%                (corner, interior etc.)
%  - nodes:      data structure containing arrays of nodes in each category
%
% HISTORY
% S.M. Giusti,         ???  2007: initial coding, based on similar source 
%                                 code written by E.A. de Souza Neto.
% S.M. Giusti,         Oct  2008: some Improvements in the code, related 
%                                 to "traction" option.
% S. Amstutz,          June 2009: modified to compatibilize with pdetool
% A.A. Novotny,        June 2009: code updating
% S.M. Giusti,         July 2009: modified to compatibilize "traction" 
%                                 option with pdetool.
% S.M. Giusti,         Feb  2011: code updating.
% A.A. Novotny,        June 2012: code updating.
% A.A. Novotny,        March 2020: code updating.
%   V. Calisti
%**************************************************************************

function [mesh,bc] = splitmesh(mesh,ndofn,bc)

% mesh and geometry  parameter
coord = (mesh.p)'; e = mesh.e;

npoin=size(coord,1); % Number of points

% multi-scale model
model = bc.model;
    
% Find maximum and minimum coordinates (coordinates of sides of
% rectangular micro-cell), dimensions and volume of micro-cell
maxcoord    = max(coord);
mincoord    = min(coord);
cell_length = maxcoord - mincoord;
bc.volRVE   = prod(cell_length);

% Find nodes on micro-cell boundary and create arrays (sets) of boundary
% and interior node numbers
% ----------------------------------------------------------------------

toler = 10^(-8)*cell_length;
maxmt = maxcoord-toler;
maxpt = maxcoord+toler;
minmt = mincoord-toler;
minpt = mincoord+toler;
node_interior = []; node_edge_left   = []; node_edge_right = [];
node_edge_top = []; node_edge_bottom = []; node_corner     = [];

isxleft   = abs( coord(:,1) - mincoord(1) ) < toler(1);
isxright  = abs( coord(:,1) - maxcoord(1) ) < toler(1);
isybottom = abs( coord(:,2) - mincoord(2) ) < toler(2);
isytop    = abs( coord(:,2) - maxcoord(2) ) < toler(2);
node_corner(1,1) = find( isybottom & isxleft  );
node_corner(2,1) = find( isybottom & isxright );
node_corner(3,1) = find(    isytop & isxright );
node_corner(4,1) = find(    isytop & isxleft  );
node_edge_left   = find(   isxleft & ~isybottom & ~isytop    );
node_edge_right  = find(  isxright & ~isybottom & ~isytop    );
node_edge_bottom = find( isybottom & ~isxleft   & ~isxright  );
node_edge_top    = find(    isytop & ~isxleft   & ~isxright  );
node_interior    = find(  ~isxleft & ~isxright  & ~isybottom  & ~isytop  );

ni = length(node_interior);
nl = length(node_edge_left);
nr = length(node_edge_right);
nb = length(node_edge_bottom);
nt = length(node_edge_top);


% Perform various checks to ensure validity of micro-cell mesh and split
% nodes into interior, "plus", "minus" and corner sets and also create
% boundary node set

%... check that given discretised rectangular cell has 4 corner nodes
corner_check = node_corner == zeros(4,1);
if any(corner_check)
    error('Wrong number of corner nodes detected in micro-cell')
end

% Check one-to-one correspondance between boundary nodes when
% periodic micro-cell boundary displacement fluctuation option is 
% selected. In this case, the sets of "plus" and "minus" boundary
% nodes must contain matching pairs of nodes on opposing faces of the
% micro-cell
if nr ~= nl
    error('Number of nodes of right and left edges of micro-cell do not coincide')
elseif nt ~= nb
    error('Number of nodes of top and bottom edges of micro-cell do not coincide')
end

% initialise array of right-left pairs
node_right_left = zeros(nr,2);
% check one-to-one correspondance and set right-left matching pairs
% array
tolery = toler(2);
[ ordered_left_y,  perm_left  ]  = sort(  coord(node_edge_left, 2 ) );
[ ordered_right_y, perm_right ]  = sort( coord(node_edge_right, 2 ) );
match  = abs( ordered_right_y - ordered_left_y ) < tolery;
if any(~match)
    error('No one-to-one correspondence between micro-cell right and left boundary nodes')
end

% initialise array of top-bottom pairs
node_right_left = [ node_edge_right(perm_right), node_edge_left(perm_left) ];
% check one-to-one correspondance and set top-bottom matching pairs
% array
node_top_bottom = zeros(nt,2);
tolerx = toler(1);
[ ordered_top_x,    perm_top    ] = sort(    coord(node_edge_top, 1 ) );
[ ordered_bottom_x, perm_bottom ] = sort( coord(node_edge_bottom, 1 ) ); 
match  = abs( ordered_top_x - ordered_bottom_x ) < tolerx;
if any(~match)
    error('No one-to-one correspondence between micro-cell top and bottom boundary nodes')
end
node_top_bottom = [ node_edge_top(perm_top), node_edge_bottom(perm_bottom) ];

[lrl,crl] = size(node_right_left);
[ltb,ctb] = size(node_top_bottom);
mesh.nodes.p(1:lrl)         = node_right_left(:,1);
mesh.nodes.p(lrl+1:lrl+ltb) = node_top_bottom(:,1);
mesh.nodes.m(1:lrl)         = node_right_left(:,2);
mesh.nodes.m(lrl+1:lrl+ltb) = node_top_bottom(:,2);

mesh.nodes.i   = node_interior';
mesh.nodes.c   = node_corner';
mesh.nodes.b   = [ mesh.nodes.p  mesh.nodes.m  mesh.nodes.c ];
mesh.nodes.all = [ mesh.nodes.i  mesh.nodes.b ];

bc.bound       = [];
bc.bound.R     = [];
bc.R           = []; 
bc.bound.Cglob = [];

% Create arrays (sets) of interior, boundary ("plus" and "minus") and
% corner nodes degrees of freedom
% -------------------------------------------------------------------

bc.dofs = struct('i',[],'p',[],'m',[],'c',[],'b',[]);

bc.dofs.i   = [mesh.nodes.i  mesh.nodes.i+npoin];
bc.dofs.p   = [mesh.nodes.p  mesh.nodes.p+npoin];
bc.dofs.m   = [mesh.nodes.m  mesh.nodes.m+npoin];
bc.dofs.c   = [mesh.nodes.c  mesh.nodes.c+npoin];
bc.dofs.b   = [bc.dofs.p  bc.dofs.m  bc.dofs.c];
bc.dofs.all = [bc.dofs.i  bc.dofs.b];
    

end
