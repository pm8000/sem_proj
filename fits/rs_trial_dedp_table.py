#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  6 17:45:43 2022

@author: pascal
"""
import CoolProp.CoolProp as CP

def get_dedp(p, rho, fluid, distance=0.001):
    #approximate derivative by linear variation in p and piecewise polynomials for rho
    #pressure should not deviate more than 1kPa from 1bar
    #pressure dependence is approximated by linear fit, with coefficents fitted depending on density
    #first based on density the coefficients for pressure dependence are evaluated
    #inputs: floats p (pressure), rho(density)
    #        string fluid(name known by CoolProp, approximation only works for water)
    #        float distance (distance bewteen the to evaluated values to compute derrivative, optional)
    #output: derivative de/dp approximated at given point
    
    if rho<0.4 or fluid!='Water':
        #not in area of interest
        return (CP.PropsSI('Umass','P',p*(1+distance),'D',rho,fluid)-CP.PropsSI('Umass','P',p*(1-distance),'D',rho,fluid))/(2*p*distance)
        
    elif rho<0.48:
        coeff=[-1.14543811e-04,  6.27731469e-05]
        linear=coeff[0]*rho+coeff[1]
        coeff=[ 53.32459343, -55.23581921,  20.09434416]
        constant=coeff[0]*rho**2+coeff[1]*rho+coeff[2]
        
    elif rho<0.52:
        coeff=[-2.72749424e-04,  1.47761711e-04, -2.02532238e-07]
        linear=coeff[0]*rho**2+coeff[1]*rho+coeff[2]
        coeff=[ 68.06858284, -69.31179626,  23.45548479]
        constant=coeff[0]*rho**2+coeff[1]*rho+coeff[2]
        
    elif rho<0.53:
        coeff=[-1.42877931e-04,  7.71797783e-05]
        linear=coeff[0]*rho+coeff[1]
        coeff=[ 87.5581132,  -89.43977943,  28.65260177]
        constant=coeff[0]*rho**2+coeff[1]*rho+coeff[2]
        
    elif rho<0.535:
        coeff=[-1.51418392e-04,  8.17006272e-05]
        linear=coeff[0]*rho+coeff[1]
        coeff=[3.89633428, 3.77913804]
        constant=coeff[0]*rho+coeff[1]
        
    elif rho<0.538:
        coeff=[-1.56954111e-04,  8.46608171e-05]
        linear=coeff[0]*rho+coeff[1]
        coeff=[4.76070368, 3.31693237]
        constant=coeff[0]*rho+coeff[1]
    
    elif rho<=0.5392:
        coeff=[-1.60305030e-04,  8.64626208e-05]
        linear=coeff[0]*rho+coeff[1]
        coeff=[5.26353409, 3.04650909]
        constant=coeff[0]*rho+coeff[1]
        
    elif rho<0.541:
        coeff=[ 2.64700e-08,  6.65084e-09, -1.31947e-08, -3.30669e-08, -5.29657e-08,
         -7.28915e-08, -9.28443e-08, -1.12824e-07, -1.32832e-07, -1.52867e-07,
         -1.72929e-07, -1.93020e-07, -2.13138e-07, -2.33285e-07, -2.53460e-07,
         -2.73664e-07]
        density=[0.5391968,  0.53931955, 0.53944229, 0.53956504, 0.53968779, 0.53981054,
         0.53993329, 0.54005604, 0.54017879, 0.54030153, 0.54042428, 0.54054703,
         0.54066978, 0.54079253, 0.54091528, 0.54103803]
        for i in range(len(coeff)-1):
            if rho>=density[i] and rho<=density[i+1]:
                interpol=(rho-density[i])/(density[i+1]-density[i])
                linear=coeff[i]+interpol*(coeff[i+1]-coeff[i])
                break
        coeff=[5.64688827, 2.83971793]
        constant=coeff[0]*rho+coeff[1]
        
    elif rho<=0.544:
        coeff=[-1.67662333e-04,  9.04393143e-05]
        linear=coeff[0]*rho+coeff[1]
        coeff=[6.28677592, 2.49330957]
        constant=coeff[0]*rho+coeff[1]
        
    elif rho<=0.548:
        coeff=[-1.76110665e-04,  9.50376715e-05]
        linear=coeff[0]*rho+coeff[1]
        coeff=[7.30336181, 1.94002977]
        constant=coeff[0]*rho+coeff[1]
        
    elif rho<=0.56:
        coeff=[-0.00209852,  0.00212287, -0.00053462]
        linear=coeff[0]*rho**2+coeff[1]*rho+coeff[2]
        coeff=[10.61072349,  0.12217741]
        constant=coeff[0]*rho+coeff[1]
        
    elif rho<=0.572:
        coeff=[-0.00493376,  0.00530396, -0.00142691]
        linear=coeff[0]*rho**2+coeff[1]*rho+coeff[2]
        coeff=[19.52752734, -4.87786207]
        constant=coeff[0]*rho+coeff[1]
        
    elif rho<=0.584:
        coeff=[-0.01212721,  0.01354645, -0.00378807]
        linear=coeff[0]*rho**2+coeff[1]*rho+coeff[2]
        coeff=[ 1270.30037682, -1428.76299685,   407.93791397]
        constant=coeff[0]*rho**2+coeff[1]*rho+coeff[2]
        
    elif rho<0.596:
        #too close to change to vapour region
        return (CP.PropsSI('Umass','P',p*(1+distance),'D',rho,fluid)-CP.PropsSI('Umass','P',p*(1-distance),'D',rho,fluid))/(2*p*distance)
    
    elif rho<0.8:
        coeff=[-3.84152726e-05,  8.02816986e-05, -6.53460273e-05]
        linear=coeff[0]*rho**2+coeff[1]*rho+coeff[2]
        coeff=[ 36.93938426, -77.19719674,  55.67757844]
        constant=coeff[0]*rho**2+coeff[1]*rho+coeff[2]
        
    elif rho<0.9:
        coeff=[ 1.77822476e-05, -3.98745557e-05]
        linear=coeff[0]*rho+coeff[1]
        coeff=[-17.09900674,  31.18476655]
        constant=coeff[0]*rho+coeff[1]
        
    elif rho<1.0:
        coeff=[ 1.42288529e-05, -3.66846206e-05]
        linear=coeff[0]*rho+coeff[1]
        coeff=[-13.68215386,  28.11740985]
        constant=coeff[0]*rho+coeff[1]
        
    elif rho<1.2:
        coeff=[-9.70246704e-06,  3.19923312e-05, -4.47661738e-05]
        linear=coeff[0]*rho**2+coeff[1]*rho+coeff[2]
        coeff=[  9.3296953,  -30.76310838,  35.8884221 ]
        constant=coeff[0]*rho**2+coeff[1]*rho+coeff[2]
        
    elif rho<2:
        coeff=[-3.34145220e-06,  1.59022181e-05, -3.45583738e-05]
        linear=coeff[0]*rho**2+coeff[1]*rho+coeff[2]
        coeff=[  3.21303625, -15.2911072,   26.07276895]
        constant=coeff[0]*rho**2+coeff[1]*rho+coeff[2]
    
    elif rho<3:
        coeff=[-8.56852862e-07,  6.38892112e-06, -2.54045854e-05]
        linear=coeff[0]*rho**2+coeff[1]*rho+coeff[2]
        coeff=[ 0.82395934, -6.14360551, 17.27102499]
        constant=coeff[0]*rho**2+coeff[1]*rho+coeff[2]  
        
    elif rho<4:
        coeff=[-3.06296424e-07,  3.20373471e-06, -2.07844436e-05]
        linear=coeff[0]*rho**2+coeff[1]*rho+coeff[2]
        coeff=[ 0.29456445, -3.08088696, 12.82858699]
        constant=coeff[0]*rho**2+coeff[1]*rho+coeff[2]
        
    elif rho<5.5:
        coeff=[-1.23267213e-07,  1.74809380e-06, -1.78818776e-05]
        linear=coeff[0]*rho**2+coeff[1]*rho+coeff[2]
        coeff=[ 0.11852644, -1.68089659, 10.03708923]
        constant=coeff[0]*rho**2+coeff[1]*rho+coeff[2]

    elif rho<7:
        coeff=[-5.34639710e-08,  9.99408857e-07, -1.58700306e-05]
        linear=coeff[0]*rho**2+coeff[1]*rho+coeff[2]
        coeff=[ 0.05138629, -0.9607178,   8.10168616]
        constant=coeff[0]*rho**2+coeff[1]*rho+coeff[2]
        
    elif rho<9:
        coeff=[ 2.02783368e-07, -1.28944157e-05]
        linear=coeff[0]*rho+coeff[1]
        coeff=[ 0.02450927, -0.5866021,   6.79729364]
        constant=coeff[0]*rho**2+coeff[1]*rho+coeff[2]
        
    elif rho<11:
        coeff=[ 1.29364935e-07, -1.22411267e-05]
        linear=coeff[0]*rho+coeff[1]
        coeff=[-0.12439249,  4.61309873]
        constant=coeff[0]*rho+coeff[1]

    elif rho<13.5:
        coeff=[ 8.61448079e-08, -1.17666749e-05]
        linear=coeff[0]*rho+coeff[1]
        coeff=[-0.08283097,  4.15685108]
        constant=coeff[0]*rho+coeff[1]
        
    elif rho<16.5:
        coeff=[ 5.73595264e-08, -1.13791185e-05]
        linear=coeff[0]*rho+coeff[1]
        coeff=[-0.05515463,  3.78421531]
        constant=coeff[0]*rho+coeff[1]
        
    elif rho<21:
        coeff=[ 3.68178644e-08, -1.10385495e-05]
        linear=coeff[0]*rho+coeff[1]
        coeff=[-0.03540437,  3.45677106]
        constant=coeff[0]*rho+coeff[1]
           
    elif rho<27:
        coeff=[ 2.24906153e-08, -1.07379514e-05]
        linear=coeff[0]*rho+coeff[1]
        coeff=[-0.02162648,  3.1677025 ]
        constant=coeff[0]*rho+coeff[1]
    
    elif rho<35:
        coeff=[ 1.34797222e-08, -1.04949382e-05]
        linear=coeff[0]*rho+coeff[1]
        coeff=[-0.0129618,   2.93402611]
        constant=coeff[0]*rho+coeff[1]
        
    elif rho<47:
        coeff=[ 7.72705727e-09, -1.02929550e-05]
        linear=coeff[0]*rho+coeff[1]
        coeff=[-0.00743056,  2.73982018]
        constant=coeff[0]*rho+coeff[1]
        
    elif rho<61:
        coeff=[ 4.44269539e-09, -1.01398620e-05]
        linear=coeff[0]*rho+coeff[1]
        coeff=[-0.00427204,  2.59259477]
        constant=coeff[0]*rho+coeff[1]
        
    elif rho<78:
        coeff=[ 2.67918117e-09, -1.00328017e-05]
        linear=coeff[0]*rho+coeff[1]
        coeff=[-0.00257624,  2.48964537]
        constant=coeff[0]*rho+coeff[1]
        
    elif rho<100:
        coeff=[ 1.63403841e-09, -9.95142787e-06]
        linear=coeff[0]*rho+coeff[1]
        coeff=[-1.57124167e-03,  2.41139748e+00]
        constant=coeff[0]*rho+coeff[1]
        
    elif rho<150:
        coeff=[ 8.41114241e-10, -9.87032406e-06]
        linear=coeff[0]*rho+coeff[1]
        coeff=[-8.08802417e-04,  2.33341181e+00]
        constant=coeff[0]*rho+coeff[1]
        
    elif rho<225:
        coeff=[ 3.74005352e-10, -9.80065587e-06]
        linear=coeff[0]*rho+coeff[1]
        coeff=[-3.59633219e-04,  2.26641931e+00]
        constant=coeff[0]*rho+coeff[1]
        
    elif rho<350:
        coeff=[ 1.59759097e-10, -9.75246946e-06]
        linear=coeff[0]*rho+coeff[1]
        coeff=[-1.53612599e-04,  2.22008249e+00]
        constant=coeff[0]*rho+coeff[1]
        
    elif rho<500:
        coeff=[ 7.23450996e-11, -9.72245404e-06]
        linear=coeff[0]*rho+coeff[1]
        coeff=[-6.95613776e-05,  2.19122122e+00]
        constant=coeff[0]*rho+coeff[1]
        
    elif rho<680:
        coeff=[ 3.73603301e-11, -9.70516602e-06]
        linear=coeff[0]*rho+coeff[1]
        coeff=[-3.59235818e-05,  2.17459827e+00]
        constant=coeff[0]*rho+coeff[1]
    
    elif rho<=958.5:
        linear=-9.67702921e-06
        constant=2.14754333
        
    elif rho<959.9:
        #too close to change to transition to pure liquid
        return (CP.PropsSI('Umass','P',p*(1+distance),'D',rho,fluid)-CP.PropsSI('Umass','P',p*(1-distance),'D',rho,fluid))/(2*p*distance)
    
    elif rho<965:
        coeff=[-4.91965464e-13,  4.59919164e-10]
        linear=coeff[0]*rho+coeff[1]
        coeff=[ 2.07995763e-05, -1.74781144e-02]
        constant=coeff[0]*rho+coeff[1]
        
    elif rho<970:
        coeff=[-6.70854264e-13,  6.32355670e-10]
        linear=coeff[0]*rho+coeff[1]
        coeff=[ 2.68136394e-05, -2.32821660e-02]
        constant=coeff[0]*rho+coeff[1]
        
    elif rho<975:
        coeff=[-1.09102932e-12,  1.03988667e-09]
        linear=coeff[0]*rho+coeff[1]
        coeff=[ 3.56703788e-05, -3.18733124e-02]
        constant=coeff[0]*rho+coeff[1]
        
    elif rho<980:
        coeff=[-1.82179578e-12,  1.75244417e-09]
        linear=coeff[0]*rho+coeff[1]
        coeff=[ 4.99631274e-05, -4.58094223e-02]
        constant=coeff[0]*rho+coeff[1]
        
    elif rho<985:
        coeff=[-3.46316830e-12,  3.36142364e-09]
        linear=coeff[0]*rho+coeff[1]
        coeff=[ 3.44790250e-06, -6.69925577e-03,  3.25706896e+00]
        constant=coeff[0]*rho**2+coeff[1]*rho+coeff[2]
        
    elif rho<990:
        coeff=[-8.17373162e-12,  8.00300800e-09]
        linear=coeff[0]*rho+coeff[1]
        coeff=[ 8.16851626e-06, -1.60026096e-02,  7.84081864e+00]
        constant=coeff[0]*rho**2+coeff[1]*rho+coeff[2]
        
    elif rho<993.5:
        coeff=[-3.03943521e-11,  3.00128991e-08]
        linear=coeff[0]*rho+coeff[1]
        coeff=[ 2.30610999e-05, -4.54931159e-02,  2.24402057e+01]
        constant=coeff[0]*rho**2+coeff[1]*rho+coeff[2]
        
    elif rho<995.8:
        coeff=[-6.63933816e-11,  6.57933348e-08]
        linear=coeff[0]*rho+coeff[1]
        coeff=[ 6.77301853e-05, -1.34250156e-01,  6.65300589e+01]
        constant=coeff[0]*rho**2+coeff[1]*rho+coeff[2]
        
    elif rho<997:
        coeff=[-6.30441767e-11,  1.25449068e-07, -6.24068050e-05]
        linear=coeff[0]*rho**2+coeff[1]*rho+coeff[2]
        coeff=[ 1.85619690e-04, -3.69015069e-01,  1.83407657e+02]
        constant=coeff[0]*rho**2+coeff[1]*rho+coeff[2]
        
    elif rho<998.24:
        coeff=[-2.51991238e-10,  5.02245013e-07, -2.50257287e-04]
        linear=coeff[0]*rho**2+coeff[1]*rho+coeff[2]
        coeff=[ 5.07166046e-04, -1.01022070e+00,  5.03069708e+02]
        constant=coeff[0]*rho**2+coeff[1]*rho+coeff[2]
        
    else:
        #not in area of interest
        return (CP.PropsSI('Umass','P',p*(1+distance),'D',rho,fluid)-CP.PropsSI('Umass','P',p*(1-distance),'D',rho,fluid))/(2*p*distance)
    
    return linear*p+constant