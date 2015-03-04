# Notes #


## Comments ##
1. `minseg_project.m`
    I tried to keep everything in one file, `minseg_project.m`, but there are two different deliverables - 
        * The report
        * The controller
    So, it seems easier to have the report generator be `minseg_project.m`
        * we don't need to use everything, but it can generate pieces we want to insert into the report
    And the controller can be run by, `run_MinSeg_Controller.m`
2. `run_MinSeg_Controller.m`
    * This can load parameters and state-space matrices into the workspace before calling the simulink model.
    * The only caveat would be to remember, before the final report to place the final parameters, etc. back into the report
3. `MinSeg_Controller.slx`
    * I'm not sure if this is correct... are the other state variables not needed to balance? 


## Good References ##
[Inverted Pendulum: State-Space Methods for Controller Design](http://ctms.engin.umich.edu/CTMS/index.php?example=InvertedPendulum&section=ControlStateSpace)

[Inverted Pendulum: System Modeling](http://ctms.engin.umich.edu/CTMS/index.php?example=InvertedPendulum&section=SystemModeling)

[Linear-Quadratic-Gaussian Control](http://www.mathworks.com/help/control/multi-input-multi-output-control-design.html)

[Design an LQG Regulator](http://www.mathworks.com/help/control/getstart/functions-for-compensator-design.html#f2-1031766)
