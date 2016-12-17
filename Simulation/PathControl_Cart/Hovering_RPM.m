function Hovering_RPM(block)
setup(block);
end

function setup(block)
 
  % Register the number of ports.
  block.NumInputPorts  = 1;
  block.NumOutputPorts = 1;
  
  % Set up the port properties to be inherited or dynamic.
  block.InputPort(1).Dimensions         = 1;
  block.InputPort(1).DirectFeedthrough = false;
  block.InputPort(1).SamplingMode      = 'Sample';
 
  %------
  block.OutputPort(1).Dimensions       = 1;
  block.OutputPort(1).SamplingMode     = 'Sample';

  block.NumDialogPrms = 1;
  block.NumContStates = 1;
  
  block.SampleTimes = [0 0];
  block.SetAccelRunOnTLC(false);
  block.SimStateCompliance = 'DefaultSimState';
  block.RegBlockMethod('CheckParameters', @CheckPrms);
  block.RegBlockMethod('Outputs', @Outputs);
end
 
function CheckPrms(block)
     drone   = block.DialogPrm(1).Data;
end

function Outputs(block)
    drone = block.DialogPrm(1).Data;
    Throttle_input = block.InputPort(1).Data;
    % Calculate Thrust
    rho = 1.25; % Air density
    rr = 0.765; % Effective air flow
    E = drone.E;
    alpha = drone.Alpha;
    beta = drone.Beta;
    D = drone.D;
    RPM_ratio = drone.b/100;
    Thrust_hovering = drone.mass*drone.g;
    
    RPM_hovering = 1000 * ((Thrust_hovering / (2*E))^3 / (rr*(pi/2)*D^2*rho*alpha^2))^(1/(2*beta));
    
    block.OutputPort(1).Data = Throttle_input + RPM_hovering/RPM_ratio;
end