

S = load("srimspline.npy")
N = load("nonion_max.npy")
#pul_corr = pandas.read_hdf("data/New_calibration.hdf",key="pulser_correction").set_index(["det","strip"])

"""
    PHD is the function that calculates the correction to E_det to get the real energy.
    it contains the follwing corrections:
    1. deadlayer
    2. non ionizing
    3. angular depandancy
    4. pulse-height-defect
    """

def phdT(Ep=3182.690,fstrip=None,bstrip=None,term="all",dl=100,linear=False,det=1):
    k = 2.8e-4*0.001 # nm/e-h pair
    eps = 3.67e-3 # keV
    from scipy.integrate import quad
    from scipy.interpolate import splev
    
    
    # term 1: deadlayer
    dl = dl*0.001 # um
    nsteps = 5
    dl_step= dl/nsteps
    Etmp = Ep
    
    if fstrip or bstrip:
        bstrip = 16.5 if bstrip==None else bstrip
        fstrip = 16.5 if fstrip==None else fstrip
        det_center_x = -1.08*25.4
        det_center_y = 0.2*25.4 # plus minus of this
        det_center_z = 4.155*25.4
        
        theta1 = numpy.arctan2(numpy.sqrt((det_center_x + (bstrip-16.5)*2)*(det_center_x + (bstrip-16.5)*2)+\
                                          (det_center_y + (fstrip-16.5)*2)*(det_center_y + (fstrip-16.5)*2))
                               ,det_center_z)
        theta2 = numpy.arctan2(numpy.sqrt((det_center_x + (bstrip-16.5)*2)*(det_center_x + (bstrip-16.5)*2)+\
                                                                 (-det_center_y + (fstrip-16.5)*2)*(-det_center_y + (fstrip-16.5)*2))
                                                      ,det_center_z)
                               
        dl_step = dl/numpy.cos(0.5*(theta1+theta2))/nsteps
    else:
        fstrip=16
        bstrip=16
    
    for i in range(nsteps):
        Etmp = Etmp - splev(Etmp,S)*dl_step
    
    dEwin = (Ep - Etmp)
    if linear==True:
        dEwin = (Ep - Etmp)
    # term 2: non ionizing
    dEn = splev(Etmp,N)
    if linear==True:
        dEn = splev(Etmp,N)
    
    
    # term 3: non linearity
    dEint = quad(lambda E: 1/(1-(k/eps)*splev(E,S)),0,Ep-dEwin-dEn)[0]
    if linear==True:
        dEint = -k/eps*quad(lambda E: splev(E,S),0,Ep-dEwin-dEn)[0]

    # term 4: adc nonlinearity
        dEadc = 0.0#poly1d(pul_corr.loc[(det,fstrip),["c","b","a"]].values)(Ep)

    if term=="all":
        return dEwin


#return Ep - (dEwin + dEn + dEint)
    #elif term==1:
#return dEwin
#elif term==2:
#return dEn
#   elif term==3:
#       return dEint
#elif term==4:
#   return dEadc
#   elif term=="1and2":
#       return dEwin + dEn
#elif term=="2and3":
#   return dEn + dEint
#   elif term=="1and2and3":
#       return dEn + dEint + dEadc

vphdT = vectorize(phdT,excluded=["term","dl","linear","det"])
