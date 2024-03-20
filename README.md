C++ implementation of the 4-point 3-view algorithm for metric reconstruction from

@article{nister2006four,<br />
&nbsp;&nbsp;&nbsp; title={Four points in two or three calibrated views: Theory and practice},<br />
&nbsp;&nbsp;&nbsp; author={Nist{\'e}r, David and Schaffalitzky, Frederik},<br />
&nbsp;&nbsp;&nbsp; journal={International Journal of Computer Vision},<br />
&nbsp;&nbsp;&nbsp; volume={67},<br />
&nbsp;&nbsp;&nbsp; number={2},<br />
&nbsp;&nbsp;&nbsp; pages={211--231},<br />
&nbsp;&nbsp;&nbsp; year={2006},<br />
&nbsp;&nbsp;&nbsp; publisher={Springer}<br />
}


Experiments were performed on the system with Intel(R) Core(TM) i5-1155G7 @ 2.5 GHz.

Notes:
- number of trials is 10^5
- optimization flag is -o2
- noise&outlier free images of synthetic scene
- error is computed as 0.5*log10(||P2 - P2gt||^2 + ||P3 - P3gt||^2)
- MAXIT is the maximum number of iterations for the golden section method
- MAXLM is the maximum number of local minima of the cost function
- cost function evaluation and local minima polishing are parallelized, # of threads is 8

1) Setup #1 (MAXLM = 5, MAXIT = 10)

	- Average runtime: 0.03 ms
	- Median error: -3.02
	- Mean error: -2.52
	- Fails (error > -2): 27.2%

2) Setup #2 (MAXLM = 10, MAXIT = 20)

	- Average runtime: 0.05 ms
	- Median error: -5.49
	- Mean error: -4.86
	- Fails (error > -2): 13.5%

3) Setup #3 (MAXLM = 20, MAXIT = 40)

	- Average runtime: 0.10 ms
	- Median error: -8.81
	- Mean error: -8.06
	- Fails (error > -2): 7.0%

4) Setup #4 (MAXLM = 50, MAXIT = 30)

	- Average runtime: 0.10 ms
	- Median error: -8.05
	- Mean error: -7.70
	- Fails (error > -2): 2.9%
