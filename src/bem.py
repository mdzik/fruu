# -*- coding: utf-8 -*-
"""
Created on Wed Jan 23 20:56:48 2013

@author: michal

Blade element method code, using QPROP Formulation, Mark Drela, MIT Aero & Astro

"""
import math

import scipy.optimize
import scipy.interpolate

import numpy

import matplotlib.pyplot as plt


class ellData:
    pass

class PowerSeeker:        
    def __call__(self, rpm, *args):
        mbo = args[0]
        pwr = args[1]            
        mbo.Omega = rpm
        mbo.performSingleAnalyse(plot=False)
        return mbo.P - pwr

class EtaSeeker:        
    def __call__(self, x, *args):

        mbo = args[0]
        cl = args[1]
        eta = args[2]

        mbo.work_element.c = x[0]
        mbo.work_element.beta = x[1]
        
        mbo.solveWE()

        return [ cl - mbo.work_element.cl, eta - mbo.work_element.eta  ]

class DesignPowerSeeker:        
    
    cl_max = 0    
    cl_max_r = 0

    
    def __call__(self, x, *args):

        eta = x
        mbo = args[0]
        Pspec = args[1]

        mbo.P = 0
        mbo.eta_avg = 0
        
        for temp in mbo.blade_elements:
            mbo.work_element = temp
            c = EtaSeeker()                        
            beta0 = math.atan(  ( mbo.V * ( ( 2. * math.pi ) / mbo.Omega ))  / (2. * math.pi * temp.r) ) 

            if( temp.rR < self.cl_max_r ):
                cls = self.cl_max
            else:
                cls = self.cl_max * 0.5 * ( 1. + ( 1.  - temp.rR ) / ( 1.0 - self.cl_max_r ))             
            

            scipy.optimize.fsolve(c , x0=numpy.array([ mbo.R / 5., beta0 ]), args=( mbo, cls, eta ) )   
            mbo.P += temp.dQ * mbo.Omega
            mbo.eta_avg += temp.eta / mbo.element_count
            
        return mbo.P - Pspec
        
class bemMain:


    
    """ Main class: stores data, computes resid and evaluates Newton interation process """

    element_count = 10
    blade_elements = list()

    def performDesign(self, **kargs):
        
        P = 0
        cl_r = list()

        
        for key in kargs:
            if key == 'V': self.V = kargs[key]
            if key == 'Omega': self.Omega = kargs[key] * math.pi / 30.       
            if key == 'P': P = kargs[key]
            if key == 'R': self.R = kargs[key]            
            if key == 'cl_max': cl_max = kargs[key]
            if key == 'cl_max_r': cl_max_r = kargs[key]
            

        beta_r = list()
        c_r = list()
        
        cl_r = list()

        c = DesignPowerSeeker()
        c.cl_max = cl_max
        c.cl_max_r = cl_max_r

        scipy.optimize.newton( c , x0=0.5, args=( self, P ) )
        
        for temp in self.blade_elements:       
            c_r.append( temp.c )    
            beta_r.append( temp.beta ) 
            cl_r.append( temp.cl )
            self.P += temp.dQ * self.Omega   

        plt.figure(3)    
        plt.subplot(311, title='Beta')
        plt.plot(self.r_r, beta_r)
        plt.subplot(312, title='Chamber')        
        plt.plot(self.r_r, c_r)
        plt.subplot(313, title='Local Lift coef.')        
        plt.plot(self.r_r, cl_r)

     #   plt.show()  
    
    def performMultipleAnalyse(self, **kargs):
        plt.figure(2)           
        RPM1 = 0
        RPM0 = 0

        V0 = 0
        V1 = 0

        mode = 'rpm' # rpm/velocity ... /power/engine
        
        
        for key in kargs:
            if ( key == 'RPM0' and self.Omega == 0 ) : RPM0 = kargs[key] * math.pi / 30.
            if ( key == 'RPM1' and RPM1 == 0 ) : RPM1  = kargs[key] * math.pi / 30.    

            if ( key == 'V0' ) : V0 = kargs[key]
            if ( key == 'V1' ) : V1 = kargs[key]

            if ( key == 'mode' ) : mode = kargs[key]
            
            if key == 'V': self.V = kargs[key]
            if key == 'Omega': self.Omega = kargs[key] * math.pi / 30.
            
            if key == 'P': P = kargs[key]
           
        N = 50
        dRPM = ( RPM1 - RPM0 ) / N
        dV = ( V1 - V0 ) / N


        Ps = list()
        Qs = list()        
        Ts = list()
        RPMs = list()
        Vs = list()        
        Etas = list()
        
        plt.figure(2)

        if ( mode == 'rpm' ):
            self.Omega = RPM0
            for i in range(0,N):
                self.performSingleAnalyse(plot=False)
                Ps.append( self.P )            
                Ts.append( self.T )  
                RPMs.append(self.Omega * 30. / math.pi )
                self.Omega += dRPM
            
            plt.subplot(211, title='Power vs Speed')
            plt.plot(RPMs, Ps)
            plt.subplot(212, title='Thrust vs Speed')        
            plt.plot(RPMs, Ts)        
           # plt.show()
        
        if ( mode == 'velocity' ):
            self.V = V0
            for i in range(0,N):
                self.performSingleAnalyse(plot=False)
                Ps.append( self.P )            
                Ts.append( self.T )  
                Vs.append(self.V )
                Etas.append(self.eta )
                self.V += dV
           # plt.figure(2)
            plt.subplot(311, title='Power vs Speed')
            plt.plot(Vs, Ps)
            plt.subplot(312, title='Torqe vs Speed')        
            plt.plot(Vs, Ts)
            plt.subplot(313, title='Eff vs Speed')        
            plt.plot(Vs, Etas)        
          #  plt.show()        
            
        if ( mode == 'power' ):
            self.V = V0
            for i in range(0,N):
                c = PowerSeeker()
                
                #c(300., self, 100.)
                scipy.optimize.newton(c , self.Omega, args=( self, P ) ) 
                RPMs.append(self.Omega * 30. / math.pi )    
                Qs.append( self.Q )   
                Ts.append( self.T )  
                Vs.append(self.V )
                Etas.append(self.eta_avg )
                self.V += dV
           # plt.figure(2)
            plt.subplot(411, title='Torqe vs Speed')
            plt.plot(Vs, Qs)
            plt.subplot(412, title='Thrust vs Speed')        
            plt.plot(Vs, Ts)
            plt.subplot(413, title='Eff. vs Speed')        
            plt.plot(Vs, Etas)      
            plt.subplot(414, title='RPM vs Speed')        
            plt.plot(Vs, RPMs)             
                      
    

            
    def performSingleAnalyse(self, **kargs):
        showPlots = False
        for key in kargs:
            if ( key == 'R' and self.R == 0 ) : self.R = kargs[key]
            if key == 'V': self.V = kargs[key]
            if key == 'Omega': self.Omega = kargs[key] * math.pi / 30.
            if key == 'plot': showPlots = kargs[key]

      #  print 'V/nD = '  + str( self.V / (self.Omega * self.R ) )
        
        self.T_r = list()
        self.Q_r = list()
        self.alpha_r = list()
        
        self.T = 0
        self.Q = 0
        
        self.solveLocals()        
        
        self.eta_avg = 0.
        
        Mas = list()
        Res = list()
        
        for temp in self.blade_elements:            
            self.Q_r.append([ temp.dT ])
            self.T_r.append([ temp.dQ ])
            self.alpha_r.append([temp.alpha * 180. / math.pi] )
            self.T += temp.dT
            self.Q += temp.dQ
            self.eta_avg += temp.eta / self.element_count
            Res.append(temp.Re)
            Mas.append(temp.Ma)
       
        self.P = self.Q * self.Omega       
        self.eta = self.V * self.T / ( self.Omega *self. Q )
      
      #
        if (showPlots):
            plt.figure(1)   
            print 'Sprawnosc = ' + str( self.V * self.T / ( self.Omega *self. Q ) )
            print 'Moc = ' + str( self. Q  * self.Omega )
            print 'Ciag = ' + str( self.T )

            plt.subplot(311, title='Local AOA')    
            plt.plot(self.r_r, self.alpha_r)

            plt.subplot(312, title='Local Ma')    
            plt.plot(self.r_r, Mas )

            plt.subplot(313, title='Local Re')    
            plt.plot(self.r_r, Res )

         #   plt.show()
        
        
    def storeLocals( self ):        
        data = numpy.hstack( ( self.r_r, self.c_r, self.beta_r, self.T_r, self.Q_r, self.alpha_r ) )   
        numpy.savetxt('output.dat', data, fmt='%f', delimiter=' ', newline='\n')        

    def __init__ (self, **kargs):       
        self.R = 0
        self.Omega = 0
        self.V = 0        
        
        self.B = 2
        self.T = 0
        self.P = 0
        self.Q = 0
        self.eta = 0
        self.ut = 0
        self.vt = 0
        self.a = 340 #m/s
        self.ua = 0
        self.ut = 0
        self.rho = 1.225
        self.element_count = 24
        
        #heppeler
        #self.mu = 0.000014607

        #drela
        self.mu = 0.17800E-04
        
        
  
        
        #drela cl/cd constants, qprop defaults
        self.cl0 = 0.5
        self.cl_a = 5.8
        self.CLmax = 1.2
        self.CLmin = -0.3

        self.clcd0 = 0.5
        self.cd0 = 0.028

        self.cd2u = 0.05        
        self.cd2l = 0.02
        
        self.Re_ref = 70000
        self.Re_ecp = -0.7
        
        self.beta_r = list()
        self.c_r = list()
        self.r_r = list()
        
    def initAir(self, T):
        R = 287.058
        p = 102400.
        lam = 1.512041288 * (10. ** -6) 
        C = 120.
        self.mu = ( lam * T**(3./2.) ) / ( T + C )
        self.rho = p / R / T
        print "Lepkosc = " + str(self.mu)
        print "Gestosc = " + str(self.rho)
        
    def loadFromFile( self ):        
        filedata = numpy.loadtxt( "input.dat", skiprows=2, delimiter="\t", usecols=range(0,7) )

        self.element_count = filedata.shape[0]
        

        beta_f = numpy.compress( [False,False,True,False,False,False,False,False], filedata, axis=1 )   
        c_f = numpy.compress( [False,False,False,False,False,True,False,False], filedata, axis=1 )  
        r_f = numpy.compress( [False,False,False,False,True,False,False,False], filedata, axis=1 )        
        
        beta_f = beta_f.reshape( beta_f.shape[0] )    
        c_f = c_f.reshape( c_f.shape[0] )    
        r_f = r_f.reshape( r_f.shape[0] )    
        
        beta_int = scipy.interpolate.interp1d(r_f, beta_f)#, kind='cubic')  
        c_int = scipy.interpolate.interp1d(r_f, c_f)#, kind='cubic')       
        
        self.R0 = r_f[0]
        self.R = r_f[r_f.shape[0] - 1]
        
        i = 0
        
        self.dr = (self.R-self.R0) / self.element_count
        
        scale = 1.
        # allocating memory for blede element data, data_row
        for i in range(0,self.element_count + 1 ):
            
            #flow-related variables, see ref, eq 17-32
            temp = ellData()

            temp.id = i
            temp.r = i * self.dr + self.R0
            
            self.r_r.append( [temp.r] )

            temp.rR = temp.r / (1.015 * self.R)
            
            #twist
            temp.beta = ( beta_int( temp.r ) ) *  math.pi / 180.

            self.beta_r.append( [temp.beta] )            
            
            #chord
            temp.c =  ( c_int( temp.r ) ) * scale
            self.c_r.append( [temp.c] )      
            #raddi


            temp.Ua = 0
            temp.Ut = 0
            temp.U = 0
            
            temp.psi = 0
            temp.Wa = 0
            temp.Wt = 0
            temp.va = 0
            temp.vt = 0
            temp.alpha = 0
            temp.W = 0
            temp.Re = 0
            temp.Ma = 0
            temp.lambda_w = 0
            temp.f = 0
            temp.F = 0

            temp.Gamma = 0
            temp.Residual = 0
            
            temp.dT = 0
            temp.dQ = 0
            temp.dD = 0
            temp.dL = 0
            
            self.blade_elements.append(temp)
            temp = 0

            

        
        foil_filedata = numpy.loadtxt( "af.dat", skiprows=6, delimiter="\t", usecols=[0,1,2] )

        alphas = numpy.compress( [True,False,False], foil_filedata, axis=1 )
        cls = numpy.compress( [False,True,False], foil_filedata, axis=1 )         
        cds = numpy.compress( [False,False,True], foil_filedata, axis=1 )   
        
        
        alphas = alphas.reshape( alphas.shape[0] )
        cls = cls.reshape( cls.shape[0] )
        cds = cds.reshape( cds.shape[0] )        

        self.cl_int = scipy.interpolate.interp1d(alphas, cls)        
        
        self.cd_int = scipy.interpolate.interp1d(alphas, cds)
#        plt.subplot(211)  
        #plt.plot(r_r, c_r )
#        
        #plt.subplot(212)  
      #  plt.plot( cls, cds )        
        #print r_r
        #plt.show()



    def computeResidual( self, psi  ):


        temp = self.work_element
        
        temp.psi = psi

        temp.Wa = 0.5 * temp.Ua + 0.5 * temp.U * math.sin( temp.psi )
        temp.Wt = 0.5 * temp.Ut + 0.5 * temp.U * math.cos( temp.psi )
        temp.va = temp.Wa - temp.Ua
        temp.vt = temp.Ut - temp.Wt
        temp.alpha = temp.beta - math.atan( temp.Wa / temp.Wt )
        temp.W = math.sqrt( temp.Wa * temp.Wa  + temp.Wt * temp.Wt )
        temp.Re = self.rho * temp.W * temp.c / self.mu
        temp.Ma = temp.W / self.a
        temp.lambda_w = ( temp.rR ) * ( temp.Wa / temp.Wt )
        
        tf = ( 0.5 * self.B * ( 1.0 - temp.rR ) / temp.lambda_w)
            
        temp.f = min( tf, 20. )    
        
        temp.F = 2.0 * math.acos( math.exp( - temp.f ) ) / math.pi
        temp.Gamma = 4.0 * temp.vt * math.pi * temp.r * temp.F * \
            math.sqrt( 1. + (4.0*temp.lambda_w*self.R  / ( math.pi * self.B * temp.r ) )**2 ) / self.B 

        temp.cl = self.getCl( temp.alpha, temp.Re, temp.Ma, temp.r )      
        
        temp.eta = self.V * temp.Wt / ( self.Omega * temp.r * temp.Wa )
        temp.dL = self.B * 0.5 * self.rho * temp.W ** 2. * self.getCl( temp.alpha, temp.Re, temp.Ma, temp.r ) * temp.c * self.dr
        temp.dD = self.B * 0.5 * self.rho * temp.W ** 2. * self.getCd( temp.alpha, temp.Re, temp.Ma, temp.r ) * temp.c * self.dr
        temp.dT = temp.dL * ( temp.Wt / temp.W ) - temp.dD * ( temp.Wa / temp.W )
        temp.dQ = ( temp.dL * ( temp.Wa / temp.W ) + temp.dD * ( temp.Wt / temp.W ) ) * temp.r
        
        temp.Residual = temp.Gamma - \
            0.5 * temp.W * temp.c * temp.cl
       
        return  temp.Residual
        
#        pme.append(temp.Gamma)

#        if k != 0 : print k / self.element_count
        
        #print residual.mean()        
        #plt.plot( residual )        

        #plt.show()
        
    def __call__( self, x, *args ):
        return self.computeResidual( x )

    
    
    def probe( self ):

        
        self.V = 15.
        self.Omega = 15000. * math.pi / 30.
        
        
        for temp in self.blade_elements:
            self.work_element = temp            
        
            temp.Ua = self.V + self.ua
            temp.Ut = self.Omega * temp.r - self.ut
            temp.U = math.sqrt( temp.Ua * temp.Ua  + temp.Ut * temp.Ut )
            
            t = list()
            x = list()
            for q in range(10,40):            
                #self.computeResidual( q / 100. )
                print temp.beta
                temp.beta = q / 20.                
                self.solveWE()
                t.append(temp.eta)                                 
                x.append(q / 10.)     
            plt.plot(x,t)
            plt.show()    
    
    def solveLocals( self ):
        for temp in self.blade_elements:
            self.work_element = temp            
            self.solveWE()
            
    def solveWE( self ):

            temp = self.work_element                   

            temp.Ua = self.V + self.ua
            temp.Ut = self.Omega * temp.r - self.ut
            temp.U = math.sqrt( temp.Ua * temp.Ua  + temp.Ut * temp.Ut )

            psi1 = math.atan2( temp.Ua, temp.Ut )
            psi2 = temp.beta + self.cl0 / self.cl_a            
            scipy.optimize.newton( self, max( psi1, psi2 ) )
        
            
    def getCl( self, alpha, Re, Ma, r ):
        beta = math.sqrt( 1. - Ma**2 )
        cl = ( self.cl0 + self.cl_a * alpha ) / beta
        if cl > self.CLmax: return self.CLmax * math.cos(alpha-self.cl0/self.cl_a)
        if cl < self.CLmin: return self.CLmin * math.cos(alpha-self.cl0/self.cl_a)  
        return cl
        #return self.cl_int( alpha * 180. / math.pi )

    def getCd( self, alpha, Re, Ma, r ):
        #return self.cd_int( alpha * 180. / math.pi )
        cl = self.getCl( alpha, Re, Ma, r )
        cd2 = self.cd2l
        if cl > self.clcd0 : cd2 = self.cd2u
        cd = ( self.cd0 + cd2 * ( cl - self.clcd0 )**2 ) * ( Re / self.Re_ref ) ** self.Re_ecp
        return cd

def main(): 
      #  plt.subplot(211)
        mB = bemMain()
        mB.loadFromFile()
#        print mB.performSingleAnalyse( R=0.0751, V=15., Omega=15000. )
        mB.initAir(271.18 + 20.)      
      #  mB.probe()
              
       #print mB.performMultipleAnalyse( R=0.0751, V=15, Omega=15000, V0=.5, V1=20, mode='power', RPM1=15000., RPM0=5000, P=30. )
        print mB.performDesign( R=0.0751, V=35., Omega=25000., cl_max=0.2, cl_max_r=0.2, P=600. )
        
        print mB.performSingleAnalyse( R=0.0751, V=35., Omega=25000., plot=True )
        
        print mB.performMultipleAnalyse( R=0.0751, V=15, Omega=25000, V0=.5, V1=60, mode='power', RPM1=15000., RPM0=5000, P=400. )
        
        plt.show()  


if __name__ == "__main__":
    main()