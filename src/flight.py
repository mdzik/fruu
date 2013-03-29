# -*- coding: utf-8 -*-
"""
Created on Sat Mar  9 12:17:49 2013

@author: michal
"""

import bem
import matplotlib.pyplot as plt
import numpy
import scipy.optimize
import scipy.interpolate



class RPMseek:

    def __init__(self, bem, engine ):
        self.bem = bem
        self.engine = engine
        
    def __call__(self, rpm ):
        self.bem.performSingleAnalyse( Omega = rpm, plot = False )
        return self.bem.T - self.engine.getTorque(rpm)

class EngineModel:
    
    def getTorque(self, rpm):
        return 0.3 * ( -((rpm-20000.)/20000)**2 + 1. )
        
    def plot(self):
        rps = numpy.arange(10., 40000.)
        power = rps * self.getTorque(rps)
        plt.figure(5)
        plt.subplot(111, title='Power vs RPS of engine')
        plt.plot(rps, power)
        #plt.show()

class flyier:
    
    
    cd = 0.0001
    rho = 0.
    mass = 0.
        
    
    time = list()
    speed = list()
    rpm = list()
    trust = list()
    drag = list()
    distance = list()
    force = list()
    
    _ode_hl = 0
    _ode_Y = list()
    _ode_F = list()
    _ode_dt = 0.75

    
    
    def getThrust(self, speed):
        self.Prop.V = self.V
        seeker = RPMseek( self.Prop, self.Engine )
        
        if (len(self.rpm) - 1 >= 0):
            
            rpms = scipy.optimize.newton( seeker , self.rpm[ len(self.rpm) - 1 ] )
        else:
            rpms = scipy.optimize.newton( seeker , self.rpm0 )
            
        self.rpm.append( rpms )
        self.drag.append(- speed*speed*self.cd)
        self.trust.append(self.Prop.T)
        self.force.append(self.Prop.T - speed * speed * self.cd )
            #    self.trust.append( self.Prop.T )
        return self.Prop.T
        
        
    def _F(self,Y):
        
        F = self.getThrust(Y[0]) - self.rho * self.cd * Y[1] * Y[1] / 2.
# F = m Y[1]
# dY[0] / dt = Y[1]
        f = numpy.array([ Y[1], F / self.mass ])
        return f
        
    def timeStep(self,F):
        Yn = self._ode_Y.pop()
        self._ode_Y.append( Yn + self._ode_dt * ( F  ) )
        
        
    def __init__(self, **kargs):
        for key in kargs:
            if key == 'cd': self.cd = kargs[key]
            if key == 'Prop': self.Prop = kargs[key]
            if key == 'V0': self.V0 = kargs[key]
            if key == 'a0': self.a0 = kargs[key]
            if key == 'Engine': self.Engine = kargs[key]
            if key == 'rpm0': self.rpm0 = kargs[key]

        self.V = self.V0
     
    def test(self):     
        for i in range(1,15):
            self.time.append(self.V)
            print self.getThrust(self.V)
            self.V = self.V + 5.
            self.speed.append(self.V)
        
        
    def solve(self):
        
        nts = 300
    #    self.rpm_int = numpy.zeros((nts,2));        
        
        for ts in range(1, nts ):
            #print ts
            V = self.V
            Drag = V * V * self.cd
            
#            if (V > self.rpm_int[ts]):            
            Thrust = self.getThrust(V)
#            else
#                intdata = 
#                t_int = scipy.interpolate.interp1d(r_f, beta_f)
#
#            self.rpm_int[ts] = [Drag, Thrust]  
#            self.rpm_int.sort(axis=1)
#            
            
            
            V = V + self._ode_dt * ( - Drag + Thrust );
            
            self.time.append(ts)
            self.speed.append(V)
           # self.drag.append(- Drag)
           # self.trust.append(Thrust)
            print V, self.V - V, ts
            self.V = V
            
            
     #   print self.rpm_int            

    
mB = bem.bemMain()
mB.loadFromFile()
mB.initAir(271.18 + 20.)      

mB.performDesign( R=0.0751, V=35., Omega=25000., cl_max=0.2, cl_max_r=0.2, P=400. )

mB.performMultipleAnalyse( R=0.0751, V=15, Omega=25000, V0=5, V1=60, mode='rpms', RPM1=35000., RPM0=5000, P=400. )


engine = EngineModel()
engine.plot()

model = flyier( Prop = mB, rpm0 = 15000., Engine = engine, V0 = 10., a0 = 1. )

model.solve()

rpm0 = numpy.arange(10000., 30000., 1000.)
plt.figure(4)
#plt.plot(rpm0, engine.getTorque(rpm0) )

plt.subplot(311)
plt.plot( model.time, model.speed )

plt.subplot(312)
plt.plot( model.time, model.rpm )

plt.subplot(313)
plt.plot( model.time, model.force )

plt.show()
#model.getThrust(50.)
    
    