# ML-Exercise

This project is based on Microsoft's Greenhouse management project, but it is simplified by my teacher Ph.D. Nguyen Tien Thinh

The problem is based on several physical equations, I had to manage to estimate the temperature in Greenhouse in the next 5 minutes, base on the Coefficient given in file 'environment.csv'. Two methods were offered to solve the problems: Explicit Euler and Runge Kutta, solved in file 'RK4-Euler.py'.

However, this exercise requires a more accurate model, so I built a deep learning model from scratch, merely NumPy module, to train and solve an ODE equation, solved in file 'ODE-solver.py'. But it's not enough to give a full solution for this exercise because it's hard to simplify the temperature equation from 20 other physical equations.
