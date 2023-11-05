% Gneral eigen problem tracking function
% INPUT:
% M(a): sqaure matrix vary with parameter a--passed by a function handle
% a0: the initial starting point of the parameter
% a_: the ending poingt of the parameter
% delta_a: stepsize of increment

%% Step 1: approximate and set up the problem
M0 = get_matrix(a0);
N = size(M0,1);
[V0,D0] =eig(M0);



% define the function getting the Jacobian matrix of the system
function Get_Jacobain()