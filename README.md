# power-flow-analysis
This MATLAB program solves for the state variables of a power system applying Newton-Raphson's algorithm

A brief explanation of this method for power flow applications can be seen here: 
https://en.wikipedia.org/wiki/Power-flow_study

A basic overview of this method
  1. Assume a voltage magnitude of 1.0pu and voltage angle of 0 degrees for all unknown variables (aka as flat start).
  2. Solve the power balance equations with the most updated voltage magnitudes and angles
  3. Linearize the system around the most recent voltage magnitudes and angles
  4. Solve for the change in voltage magnitudes and angles
  5. Replace the previous voltage magnitudes and angles values with the newly computed values
  6. Check if the desired tolerance level has been met which would indicate the solution has converged, else go back to step 2
