With find_r_candidaters.py, we find the representatives for R, as in the thesis. 

In main_calculations.py, we check whether a triple (P, Q, R) defines a (smooth) 
curve of genus 5 that is a complete intersection of three quadratic hypersurfaces.

That is a code dealing with the case when P = P6, when P6_list_of_Qs has 134 elements. 
In code, we manually change LOWER_BOUND and UPPER_BOUND after the finished job. 
Theoretically, we would set that LOWER_BOUND = 0, and that UPPER_BOUND = 134 (the length
of P6_list_of_Qs) if we would have "super-fast" computer. 

The codes for P1, P3 and P4 are very similar to this one. We only change the initial data 
contained in prepare_initial_data.py and the upper/lower bounds in the main program. 