#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  5 15:45:33 2022

@author: pascal
"""

import numpy as np
import CoolProp.CoolProp as CP

def get_dedrho(rho, fluid, distance=0.001):
    #pressure fixed at 1 bar
    #piecewise cubic approximation of derrivative
    
    if rho<0.4:
        #density not in domain (T higher than 270°C)
        return (CP.PropsSI('Umass','P',1e5,'D',rho*(1+distance),fluid)-CP.PropsSI('Umass','P',1e5,'D',rho*(1-distance),fluid))/(2*rho*distance)
    
    elif rho<0.5:
        return 9.20328450e+07*rho**3-1.53782378e+08*rho**2+9.01026774e+07*rho-1.93829389e+07
    
    elif rho<=0.5875:
        return 16920652.18149624*rho**3-42056094.17118616*rho**2+34589494.80004115*rho-10169056.44019231
    
    elif rho<0.5935:
        #change to vapour area, highly dependant on pressure
        return (CP.PropsSI('Umass','P',1e5,'D',rho*(1+distance),fluid)-CP.PropsSI('Umass','P',1e5,'D',rho*(1-distance),fluid))/(2*rho*distance)
    
    elif rho<0.7:
        return 4.40996111e+07*rho**3-1.06919346e+08*rho**2+9.20799115e+07*rho-2.97087252e+07
    
    elif rho<0.9:
        return 15457301.69543689*rho**3-46336342.88215045*rho**2+49277215.19442039*rho-19607675.54479533
    
    elif rho<1.2:
        return 3593868.42053682*rho**3-14469002.8917335*rho**2+20633440.38598556*rho-10991832.74422337
    
    elif rho<1.5:
        return 1123340.88234657*rho**3-5683458.75498856*rho**2+10205976.9937587*rho-6860541.8418124
    
    elif rho<2:
        return 311108.08836213*rho**3-2039606.13283557*rho**2+4740444.41361017*rho-4119503.69733582
    
    elif rho<2.5:
        return 87350.98701097*rho**3-736576.25464451*rho**2+2204491.03064739*rho-2469795.06304912
    
    elif rho<3.2:
        return 26912.9052097*rho**3-287418.6925272*rho**2+1089029.53601937*rho-1544024.25950675
    
    elif rho<4:
        return 8330.43928249*rho**3-112392.61698045*rho**2+538205.81802838*rho-964763.69650266
    
    elif rho<5.8:
        return 1849.28435858*rho**3-33924.55843487*rho**2+220205.20238375*rho-533387.08872016
    
    elif rho<8:
        return 329.27184985*rho**3-8509.31638247*rho**2+77903.92724421*rho-266476.39767171
    
    elif rho<11:
        return 6.65003152e+01*rho**3-2.36618415e+03*rho**2+2.98282855e+04*rho-1.40499125e+05
    
    elif rho<15:
        return 1.38290219e+01*rho**3-6.73383971e+02*rho**2+1.16189793e+04*rho-7.49236874e+04

    elif rho<22:
        return 2.41914043e+00*rho**3-1.67534213e+02*rho**2+4.10411571e+03*rho-3.75063397e+04
    
    elif rho<32:
        return 3.64402499e-01*rho**3-3.68338985e+01*rho**2+1.31729144e+03*rho-1.75784952e+04
        
    elif rho<45:
        return 6.12049775e-02*rho**3-8.82415182e+00*rho**2+4.50495349e+02*rho-8.58904790e+03
            
    elif rho<65:
        return 1.03672832e-02*rho**3-2.13479652e+00*rho**2+1.55558396e+02*rho-4.23034331e+03
                
    elif rho<90:
        return 1.84389294e-03*rho**3-5.35199470e-01*rho**2+5.50278710e+01*rho-2.11371390e+03
                    
    elif rho<130:
        return 3.23977599e-04*rho**3-1.33424783e-01*rho**2+1.94447995e+01*rho-1.05758583e+03
                        
    elif rho<185:
        return 5.35994933e-05*rho**3-3.16100138e-02*rho**2+6.59926407e+00*rho-5.14371267e+02

    elif rho<250:
        return 1.05249836e-05*rho**3-8.57505725e-03*rho**2+2.47611252e+00*rho-2.67260607e+02
    
    elif rho<350:
        return 2.12787022e-06*rho**3-2.39060206e-03*rho**2+9.51143539e-01*rho-1.41341163e+02
        
    elif rho<480:
        return 4.17740037e-07*rho**3-6.49327620e-04*rho**2+3.57603312e-01*rho-7.35919285e+01
            
    elif rho<650:
        return 8.90226940e-08*rho**3-1.88407881e-04*rho**2+1.41317540e-01*rho-3.96191281e+01
                
    elif rho<=958.6:
        return 1.56021952e-08*rho**3-4.69735040e-05*rho**2+5.00184558e-02*rho-1.98660394e+01
    
    elif rho<959:
        #change to liquid area, highly dependant on pressure
        return (CP.PropsSI('Umass','P',1e5,'D',rho*(1+distance),fluid)-CP.PropsSI('Umass','P',1e5,'D',rho*(1-distance),fluid))/(2*rho*distance)

    elif rho<980:
        return -5.31908131e-02*rho**3+1.52798371e+02*rho**2-1.46365772e+05*rho+4.67460458e+07
    
    elif rho<990:
        return -6.30156844e-01*rho**3+1.85105780e+03*rho**2-1.81261388e+06*rho+5.91696655e+08
        
    elif rho<995:
        return -7.15865661e+00*rho**3+2.12524163e+04*rho**2-2.10315946e+07*rho+6.93781296e+09
            
    elif rho<997:
        return -6.23387867e+01*rho**3+1.85975495e+05*rho**2-1.84941464e+08*rho+6.13047269e+10
    
    elif rho<998.9:
        return -770.7870983507996*rho**3+2305857.0560457795*rho**2-2299373844.4281898*rho+764303845240.5865
    
    else:
        #no longer in area of interest
        return (CP.PropsSI('Umass','P',1e5,'D',rho*(1+distance),fluid)-CP.PropsSI('Umass','P',1e5,'D',rho*(1-distance),fluid))/(2*rho*distance)