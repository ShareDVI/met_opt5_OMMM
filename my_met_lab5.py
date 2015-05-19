#!/usr/bin/python
import math

def dihotomia(function,a,b,eps=1e-8):
	iterations = 0
	left = a
	right = b
	center = (left + right)/2
	length = right - left
	while length > eps:
		an = left + length/4
		bn = right - length/4
		if function(an) < function(center):
			right = center
			center = an
		elif function(bn) < function(center):
			left = center
			center = bn
		else:
			left = an
			right = bn
		length = right - left
		iterations+=1
	return center,iterations

def gold(function,a,b,eps=1e-8):
	iterations = 0
	y = a + (3.0 - math.sqrt(5.0))*(b-a)/2
	z = a + b - y
	while b - a > eps:
		if function(y) <= function(z):
			b = z
			z = y
			y = a + b - y
		else:
			a = y
			y = z
			z = a + b - z
		iterations+=1
	return (a+b)/2,iterations
	
def oom(function,a,b,tau=0.4,eps=1e-8):
	iterations = 0
	while b - a > eps:
		y = a + tau*(b-a)
		z = a + b - y
		if function(y) <= function(z):
			b = z
		else:
			a = y
		iterations+=1
	return (a+b)/2,iterations
	

	
def func(x):
	return x**2

def exp(x):
	return math.exp(x)

if __name__ == "__main__":
	print dihotomia(exp, -10.0, 10.0),gold(exp, -10.0, 10.0),oom(exp, -10.0, 10.0)
	

	
