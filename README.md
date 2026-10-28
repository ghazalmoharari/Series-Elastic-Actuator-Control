# Series Elastic Actuator (SEA) Control System Design

## Overview
This repository contains a comprehensive two-phase control systems project focused on the modeling, time-domain analysis, and frequency-domain control of a Series Elastic Actuator (SEA)[cite: 8, 2]. SEAs are widely used in interactive robotics to increase safety and absorb shocks, but their inherent flexibility introduces unwanted vibrations and control challenges[cite: 8]. 

The system consists of a motor with inertia $J_1$ driving a robotic arm with inertia $J_2$ through a non-linear torsional spring characterized by the equation $T_{spring} = k(\theta_1 - \theta_2) + \beta(\theta_1 - \theta_2)^3$[cite: 8]. 

## Project Structure

### Phase 1: Time-Domain Modeling and Control
The first phase focuses on establishing the mathematical foundation and stabilizing the system in the time domain[cite: 7].
* **Non-linear Modeling & Linearization:** Derived the non-linear differential equations and linearized the system around the operating point $\alpha_0 = 0.1$ rad to find the equivalent stiffness[cite: 9].
* **State-Space & Transfer Function:** Extracted the state-space representation and transfer function $G(s) = \frac{\theta_1(s)}{\tau(s)}$, analyzing the physical significance of the system's anti-resonance zeros[cite: 9].
* **Control Design:** Designed a closed-loop controller with unity negative feedback to achieve[cite: 10]:
  * Absolute stability[cite: 10].
  * Overshoot $< 20\%$[cite: 10].
  * Settling time $< 2$ seconds[cite: 10].
  * Zero steady-state error for a step input[cite: 10].
* **Disturbance Rejection:** Validated the controller's robustness against a 2 Nm triangular torque disturbance applied directly to the arm[cite: 10, 11].

### Phase 2: Frequency-Domain Analysis and Compensation
The second phase shifts to frequency-domain techniques to handle resonance and improve response speed[cite: 2, 3].
* **Bode Plot Analysis:** Analyzed the linearized system's frequency response to identify resonance and anti-resonance frequencies, alongside their severe phase shifts[cite: 3].
* **Notch & Lead Compensator Design:** 
  * Designed a Notch filter ($N(s) = \frac{s^2 + 2\zeta_z\omega_0s + \omega_0^2}{s^2 + 2\zeta_p\omega_0s + \omega_0^2}$) to suppress the anti-resonance peak while maintaining overall stability[cite: 3].
  * Implemented a Lead compensator ($C_{lead}(s) = K\frac{1 + s/\omega_z}{1 + s/\omega_p}$) to improve the phase margin to at least 45 degrees and increase the bandwidth[cite: 3, 4].
* **Non-linear Validation:** Applied the final frequency-domain controller back to the original non-linear plant to verify stability, overshoot, and settling time under realistic conditions[cite: 4].

## Software & Tools
* **[Software 1, e.g., MATLAB]**: Used for calculating transfer functions, roots, and Bode plots.
* **[Software 2, e.g., Simulink]**: Used for block diagram representation and non-linear step-response simulations.

## How to Run
1. Clone this repository.
2. Open the `[Main_Script_Name.m]` in [MATLAB].
3. Run the script to generate the step responses, Bode plots, and disturbance rejection graphs for both Phase 1 and Phase 2.
