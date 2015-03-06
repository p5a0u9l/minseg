# Notes #

## Git Stuff ##
1. To see which m-files and simulink files are being tracked by git on branch `master` issue
    * `git ls-tree -r master | grep -E "(*\.m|*.slx)" `

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

## Emails ##
1. 
    From: Schulenberg, Christopher S 
    Sent: Thursday, March 05, 2015 9:09 PM
    To: Adams, Paul R
    Subject: Running the minseg untethered

    I was able to get MATLAB 2014b installed with the Arduino support package and run the model un-tethered.

    1.  Use the ‘deploy to target hardware’ button to deploy the model to Arduino:
     
    2.  Do not press the ‘play’ button
    3.  Untether the minseg
    4.  Power on the minseg

    I think when you press the play button it loads the ‘external mode’ code which relies on triggers from PC Simulink to run. If you deploy to target HW it is ready to go at the next power cycle.

    The code with no modifications gives me ~2 seconds before scooting off to oblivion. Tomorrow’s task for me is to see how modify things to 1) reflect the analytical portion of the project in order to 2) actually stabilize the minseg, and 3) clean up the god awful project directory that was unzipped from the download.

    I’d also like to include a few features like cutting the motor when the thing is laying on its side.

    Thanks,
    Chris
