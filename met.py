#!/usr/bin/python
import numpy as np
import math

def optimize_1d(f):
    alpha = 1.
    while f(alpha)>f(0):
        alpha /= 2
    return alpha

def optimize_1d_parabolas(f, a=0., b=1., eps=1e-5, feps=1e-10):
    if math.fabs(b-a)<=eps:
        return (a+b)/2
    x1=min(a,b)
    x3=max(a,b)
    c_eps=eps/(x3-x1)
    division_const=50*c_eps
    c=0.5

    #Preparation - we need to select x2
    #to be a convex down point
    #We just look at middle, and then go here and there with step=division_constant
    x2=c*x1+(1-c)*x3
    while c_eps<c<1-c_eps and (f(x2)>c*f(x1)+(1-c)*f(x2)+feps):
        if c<=0.5: c-=division_const
        else: c+=division_const
        c=1-c
        x2=c*x1+(1-c)*x3
        print "f(",x1,")=",f(x1),"; f(",x2,")=",f(x2),"; f(",x3,")=",f(x3)
    if c<=c_eps or c>=1-c_eps: #function is linear or convex up - just return min of endpoints
        if f(x1)<f(x3):
            return x1
        else:
            return x3

    #Ok, so here we got a convex down function
    i=0
    u=x2
    print "Optimizing using parabolas"
    print "i\tx1\t\t\tx2\t\t\tx3\t\t\tf(x1)\t\t\tf(x2)\t\t\tf(x3)"

    while math.fabs(x2-x1)>=feps and math.fabs(x3-x2)>=feps:
        print i, "\t",x1, "\t",x2, "\t", x3, "\t", f(x1), "\t", f(x2),"\t", f(x3)
        i=i+1

        if math.fabs((x2-x1) * (f(x2)-f(x3)) -(x2-x3) * (f(x2)-f(x1))) < eps*feps:
            #our parabola is linear...
            if f(x1)<=f(x3):
                return x1
            else:
                return x3

        u= x2 - ((x2-x1)**2 * (f(x2)-f(x3)) -(x2-x3)**2 * (f(x2)-f(x1)))/(2*((x2-x1) * (f(x2)-f(x3)) -(x2-x3) * (f(x2)-f(x1))))
        if x1<u<=x2 and f(u)-f(x1)<=feps and f(u)-f(x2)<=feps:
            x3=x2
            x2=u
        elif x2<u<x3 and f(u)-f(x2)<=feps and f(u)-f(x3)<=feps:
            x1=x2
            x2=u
        elif u<=x1 and f(u)-f(x1)<feps:
            return x1
        elif u>=x3 and f(u)<f(x3)<feps:
            return x3
        else:
            if f(x1)-f(x3)<feps:
                return x1
            else:
                return x3
    return u

class Function:
    def __init__(self, f, N):
        self.f = f
        self.N = N

    def get_gradient(self, x):
        eps=0.000001
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
            print "#",iter,"xk=",xk,"alpha=",alpha,"nullity=",nullity
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
print F.get_minimum(np.matrix([5.1, 1.1]),optimize_1d)
print F.get_minimum(np.matrix([5.1, 1.1]),optimize_1d_parabolas)

