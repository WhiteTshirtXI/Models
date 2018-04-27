clear
%  2D Stokes code with an 
%  with immersed boundary worm using implicit-time stepping

% add path for the source code
%
addpath('./src/');

% info for a restart
%
restart = false;
lastsaved = .5;

% define the domain and the grid spacing
%
Lx = 2;
Ly = 2;

% center the domain at the origin
%
xmin=-Lx/2;
ymin=-Ly/2;

% number of grid points and the grid spacing
%
K  = Lx/Ly;
Ny = 128;
Nx = K*Ny;
dx = Lx/Nx;

%fluid parameters
lam = 0;
xi = 0.5;
diffconst = 8;
%tend = 10.0;
%t0 = 0;                %<--- i think this is unnecessary
%savetime = 1;


%initialize Shat
[params,Shat,newRHS] = get_4roll_inputs(Ny,lam,xi,diffconst);


%unpack params
nu = params.nu;

% worm paramters
%
%
L  = 1.2;                      % length of the worm [mm]
Tper = 0.5;                    % period of the gait [s]
Nprefact=1.0;                  % ratio of ds to dx
N  = round( L/(Nprefact*dx) ); % number of points to discretize the worm
ds = L/(N-1);                  % length of each segment
s  = ds*(0:(N-1))';            % arclength coordinate of the worm
 
kb = 40;                        % bending stiffness
ks = 2500;                     % stretching stiffness

% time stepping information
%
dt     = 1e-3;           % time step [s]
Tend   = 7.5;  % end time [s]

Nt     = round(Tend/dt);  % number of time steps to take
saveit = round((Tper/100)/dt);  % frequency of output swimmer positions
saveall = round((Tper/50)/dt);      % frequency of output all data

% solver tolerances
%
rtol     = 1e-4;    % relative tolerance for objective function
rXtol    = 1e-12;   % relative tolerance for function values 
gmrestol = 5e-5;    % gmres tolerance for jacobian solve
  

% output locations
%
datadir    = './data';  
runname    = '2D_K_KB_40_DE_0';
fileprefix = sprintf('%s_n%03d',runname,Ny);
paramfile  = sprintf('%s/PARAMS_%s.txt',datadir,fileprefix);


%%%%%%%%%%-----END INPUT PARAMETERS-----%%%%%%%%%%


% pack up the grid parameters into a data structure
%
grid.Lx   = Lx;
grid.Ly   = Ly;
grid.xmin = xmin;
grid.ymin = ymin;
grid.Nx   = Nx;
grid.Ny   = Ny;
grid.N    = N;
grid.ds   = ds;
grid.dx   = dx;

% preallocate for speed
%
outcount_tot=round(Tend/(saveit*dt))+1;
% XTworm=zeros(N,2,outcount_tot);


% initialize swimmer body position
%
kappa0 = kappa_kicker(s,0,Tper,L);
% kappa0 = kappa_burrower(s,0,Tper,L);
X = initialize_worm(kappa0,ds);
% write parameters to file
%
fileID = fopen(paramfile,'w');
fprintf(fileID,'Lx = %f\n',Lx);
fprintf(fileID,'Ly = %f\n',Ly);
fprintf(fileID, 'Ny = %d\n',Ny);
fprintf(fileID,'dt = %0.8f\n',dt);
fprintf(fileID,'N (worm points) = %f\n',N);
fprintf(fileID, ' bending stiffness, kb = %4.4f\n',kb);
fprintf(fileID, ' stretching stiffness, ks =  %4.4f\n',ks);
fprintf(fileID, ' End time = %f\n',Tend);
fclose(fileID);

% initialize output counter and time -- overwritten if restart
%
outcount=1;
t = 0.0;

% if this is a restart, reload from the appropriate data file
%
if(restart)
    t=lastsaved;
    outcount=round(t*round(1/(saveit*dt))+1)
    fout1 = sprintf('%s/%s_t%f.mat',datadir,fileprefix,t);
    load(fout1);
    Nt   = round((Tend-t)/dt);
    X = XTworm(:,:,outcount);
end
   
% begin main loop in time
%
tic
endStep = Nt;
for tint=0:Nt
    fprintf('Time step %i of %i, time=%f \n',tint,Nt,t);
    % record the worm positons
    %
    if(mod(tint,saveit)==0)
        XTworm(:,:,outcount)=X;
        outcount=outcount+1;
        fprintf('  recording output number %g \n',outcount-1);
        if(lam~=0)
            twoNorm(outcount) = stress_2norm(Shat); 
            maxNorm(outcount) = stress_max_norm(Shat);
            fprintf(' computing stress norms %g \n', outcount-1);
        end
    end
    
    % save all the data 
    %
    if( mod(tint,saveall)==0 && t~=0)
        foutw = sprintf('%s/%s_t%f.mat',datadir,fileprefix,t);
        save(foutw,'U','Uw','XTworm','Shat');
        if(lam~=0)
            foutw2 = sprintf('%s/%s_2norm_t%f_to_t%f.mat',datadir,fileprefix,0,t);
            foutw3 = sprintf('%s/%s_maxnorm_t%f_to_t%f.mat',datadir,fileprefix,0,t);
            save(foutw2,'twoNorm');
            save(foutw3,'maxNorm');
        end
    end
    
     % New exit step
    if(tint >= endStep)
        break;
    end   
    
    % compute the curvature at the current time
    %
    kappa0 = kappa_kicker(s,t,Tper,L);
    % kappa0 = kappa_burrower(s,t,Tper,L);
    % set the external body forces to zero for stokes solve
    %   this coudl be nonzero for viscoelastic fluid
    %
    if(isnan(Shat(1,1,1)))  % if the stress blows up stop running
        break;
    end  
    if(lam == 0)
        fbhat = zeros(Nx,Ny,2);
    else  % Compute VE force
         fbhat = get_veforcehat(Shat,xi,grid);
    end

    % advance swimmer in time using backward euler
    %pass out newRHS Shat// pass in lam nu Shat newRhS
    [X,Uw,U,output] = IMstep_stokes_newton(X,dt,fbhat,ks,kb,kappa0,grid,rtol,rXtol,gmrestol); %do not need to pass Shat onwards
    if(lam~=0)
        Uhat = fft2(U);
        [Shat, newRHS] = update_Shat(Uhat,grid,Shat,nu,dt,lam,newRHS);
    end
    % update time
    %
    t = t+dt;
    %toc
end

comp_time = toc;

fprintf('total computation time  = %g \n',comp_time);
