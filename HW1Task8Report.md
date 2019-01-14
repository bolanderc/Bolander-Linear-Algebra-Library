# Homework 1 Task 8: Disaster Articles Summary



This file will give a brief summary and commentary on the three disaster articles on [this site](http://www-users.math.umn.edu/~arnold//disasters/). 

## The Patriot Missile Failure

The American Patriot Missile battery in Dharan, Saudi Arabia failed to track and intercept an incoming Iraqi missile in February of 1991 due to a numerical round-off error. Specifically, the missile system kept track of time using an integer value that represented the number of 1/10th seconds that had elapsed since the system was turned on. Error was introduced in the conversion from an integer to a floating-point number. In addition, the code used for the system was modified from an older code that was not implemented for the specific system. The new implementations were not consistent throughout the code, which allowed further inaccuracies to propagate.

## The Explosion of the Ariane 5

This article covers the explosion of a European Space Agency rocket, the Ariane 5, that self-destructed 40 s into its first voyage due to another conversion error. A software bug in the rocket's inertial reference system caused the self-destruct protocol to activate. The bug was in trying to convert a 64 bit floating point number representing the horizontal velocity of the rocket into a 16 bit signed integer value. The system and its backup crashed and sent false data through the pipeline, causing the on board computer to make an unnecessary maneuver to correct. This maneuver threatened to rip the boosters off of the rocket, which in turn triggered the self-destruct mechanism.

##  The Sinking of the Sleipner A Offshore Platform

This article talks about the Sleipner A, an oil rig in the North Sea that was destroyed after a failure in a cell wall caused a crack and leakage that the facility could not cope with. The cell wall was cracked because of a numerical error in an approximation used by the finite element program used to design the tricell (a bracing mechanism between the cells). The shear stresses in the tricell were underestimated by 47% in the finite element approximation, which led to the failure. Additional details on the exact numerical issues with the code are not examined in the summary on the website or in the papers posted.