#!/usr/bin/python
import numpy as np
import math

def optimize_1d(f):
	alpha = 1.
	while f(alpha)>f(0):
		alpha /= 2
	return alpha

class Function:
	def __init__(self, f, N):
		self.f = f;
		self.N = N;

	def get_gradient(self, x):
		eps=0.000001;
		res = np.zeros((1,self.N))
		for i in range(self.N):
			h = np.zeros((1,self.N))
			h[0,i]=eps
			res[0, i] = (self.f(x + h) - self.f(x - h)) / 2 / eps
		return res

	def get_minimum(self, x0, min_1d=optimize_1d):
		xk = x0;
		h = -1*self.get_gradient(x0);

		nullity = 1;
		iter = 0;
		while True:
			iter+=1
			if iter > 500:
				break
			f_alpha = lambda alpha: self.f(xk+alpha*h)
			alpha = min_1d(f_alpha)
			
			xk1 = xk + alpha * h
			nullity = np.linalg.norm(alpha * h)
			if nullity < 1e-10:
				break
			
			g1 = self.get_gradient(xk)
			g2 = self.get_gradient(xk1)
			beta = np.vdot(g1,g2-g1)/(np.linalg.norm(g1)**2)
			
			xk = xk1
			h = -1 * g1 + beta * h;

		return xk;

def exp2(x):
	a = x[0,0]
	b = x[0,1]
	return math.atan(a*a+b*b+2*a+2*b)

F = Function(exp2, 2)
print F.get_minimum(np.matrix([5.1, 1.1]))
