xdel(winsid());
clear;

// ------------------------- Function Definition ------------------------- //
function Id = Diode_func(Vd)
//    // With Hard Current Limiter of 10 Ampere
//   Id = min(Is.*(exp(Vd./VT) - 1),10);
    
    // Without Hard Current Limiter
    Id = Is.*(exp(Vd./VT) - 1);
endfunction

function dX = System_func(t,X)
    // Fetch Parameter
    Vout = X(1,:);
    
    // Input Signal : Intermediate Variable
    vx = sin(2*%pi*fin*t);
    
    // Input Signal : Sinusoidal Wave
    vin = Vin*vx;
    
//    // Input Signal : Square Wave (Fast Transition))
// vin = Vin*sign(vx);
//    
//    // Input Signal : Square Wave (Slow Transition)
//    kx = 20;
//    vin = Vin*(1 - exp(-kx.*vx))./(1 + exp(-kx.*vx));
    
    // Diode Current
    Vd = vin - Vout;
    Id = Diode_func(Vd);
    
    // Rectifier Model
    Gc = Id/C - Vout/(R*C);
    dX = Gc;
endfunction

function X = RK4_func(Npts,Xo)
    // Fetch State Number
    Ns = size(Xo,1);
    
    // Initial Condition
    X = zeros(Ns,Npts);
    X(:,1) = Xo;
    tn = 0;
    
    for npts = 1:Npts-1
        // Initialization
        Xn = X(:,npts);
                
        // Solver
        Ma = System_func(tn,Xn);
        Xa_est = Xn + TsampSys*Ma/2;
        Mb = System_func(tn + TsampSys/2,Xa_est);
        Xb_est = Xn + TsampSys*Mb/2;
        Mc = System_func(tn + TsampSys/2,Xb_est);
        Xc_est = Xn + TsampSys*Mc;
        Md = System_func(tn + TsampSys,Xc_est);
        Mavg = (Ma + 2*Mb + 2*Mc + Md)/6;
        X(:,npts+1) = Xn + TsampSys*Mavg;
        tn = tn + TsampSys;
    end
endfunction

// ------------------- Definition of Simulator System ------------------- //
// Samping Period
TsampSys = 1e-6;

// Simulation Interval
Tstart = 0;
Tend = 1e-2;

// Arrange Tstart and Tend
Tsim = [Tstart Tend];

// Calculate Number of Points
N = 1 + round((Tend-Tstart)/TsampSys);

// Clear Variables
clear Tstart Tend

// Time Points
Ts = Tsim(1):TsampSys:Tsim(2);

// ------------------------- System Definition -------------------------- //
// Load Parameter
fo = 12e1;
R = 250;
C = 1/(2*%pi*fo*30e-2);

// Diode Characteristics
VT = 1/40;
Is = 1e-13;

// Input Signal
fin = 1e3;
Vin = 10;

// -------------------------- Rectification --------------------------- //
// Initital Condition and ODE Solver
Xo = 0;
Xout = RK4_func(N,Xo);

// Fetch Parameter
Vout = Xout(1,:);

// -------------------------------- Plot -------------------------------- //
// Plot Data : Output Voltage over Time
figure()
plot(Ts,Vout);
