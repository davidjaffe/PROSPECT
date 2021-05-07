#!/usr/bin/env python
'''
find the reactor 

detector center to reactor center from PRD Figures 1, 2
reactor center at (0,0,0)
detector center at x, y, z = 5970 , 5090, -1190
z length of module 1176 mm
x length of detector = 2040 mm
y length of detector = 1602.85 mm (calculated assuming 14 x 11 modules in x v y)

units: mm, MeV

20210312
'''
import os,sys,datetime
import math,numpy,collections

from scipy.optimize import minimize
import Logger
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,  AutoMinorLocator)

class findreactor():
    def __init__(self,debug,TotalExpt=10,drawToFile=False,nFitPar=4,GenEvtsPerExpt=1000000,useChi2=False,setNZ=1):
        f1 = ' {0:.1f} '
        f2 = ' {0:.1f} {1:.1f} '
        f3 = ' {0:.1f} {1:.1f} {2:.1f} '
        n3 = ' {0:0.0f} {1:0.0f} {2:0.0f}'
        now = datetime.datetime.now()
        self.now = now.strftime('%Y%m%dT%H%M%S')
        
        self.drawToFile = drawToFile
        self.figDir = 'FIGURES/findreactor/'+self.now

        if not os.path.exists(self.figDir):
            os.makedirs(self.figDir)
            print 'findreactor.__init__ create self.figDir',self.figDir
        else:
            print 'findreactor.__init__ overwrite self.figDir',self.figDir

        lf = self.figDir + '/logfile.log'
        sys.stdout = Logger.Logger(fn=lf)
        print 'findreactor.__init__ Output directed to stdout and',lf

        self.nFitPar = nFitPar
        self.GenEvtsPerExpt = GenEvtsPerExpt

        self.useChi2 = useChi2
        self.setNZ   = setNZ
        
        self.debug = debug
        self.TotalExpt = TotalExpt
        print 'findreactor.__init__ debug',debug,'TotalExpt',TotalExpt,'drawToFile',drawToFile,'nFitPar',nFitPar,'GenEvtsPerExpt',GenEvtsPerExpt,'useChi2',self.useChi2,'setNZ',self.setNZ
        self.Xreactor = [0.,0.,0.]
        self.XreactorAlt = [[0.,0.,-1200.],
                            [0.,-500., 0.],
                            [-500., 0., 0.],
                            [-5000., 0., 0.],
                            [0., -5000., 0.]
                                ]
        self.nstart = nstart = 10
        maxdist = 1000.
        self.startPoints = []
        for i in range(nstart):
            self.startPoints.append( maxdist*(numpy.random.random(3)-0.5) )
        print 'findreactor.__init__ self.startPoints',self.startPoints
        self.XreactorAlt = [] ##### SET TO ZERO LENGTH #####
        self.Xdet= [5970., 5090., -1190.]     
        self.dX = [2040., 1602.85, 1176.]
        self.nx = 14
        self.ny = 11
        self.nz = self.setNZ
        self.totalModules = self.nx * self.ny * self.nz
        self.nxyz = numpy.array( [float(self.nx), float(self.ny), float(self.nz)] )
        print 'findreactor.__init__ nxyz'+n3.format(*self.nxyz)
        if self.nz>1 : print 'findreactor.__init__ Hey, NZ is NOT unity!!!!!!!'
        self.dXmod = numpy.array( self.dX ) / self.nxyz
        print 'findreactor.__init__ dXmod'+f3.format(*self.dXmod)
        print 'findreactor.__init__ Xreactor'+f3.format(*self.Xreactor),'Xdetector'+f3.format(*self.Xdet)
        xc = numpy.array(self.Xdet)

        # detector maxima, minima coordinates
        self.maxXdet = xc + 0.5*numpy.array(self.dX)
        self.minXdet = xc - 0.5*numpy.array(self.dX)
        xmi,xma = self.minXdet,self.maxXdet
        print 'findreactor.__init__ xma'+f3.format(*xma),'xmi'+f3.format(*xmi)
        xr = numpy.array(self.Xreactor)
        # min,max baselines and center-to-center baseline
        rlo = numpy.minimum(abs(xmi-xr),abs(xma-xr))
        rhi = numpy.maximum(abs(xmi-xr),abs(xma-xr))
        print 'findreactor.__init__ rlo'+f3.format(*rlo),'rhi'+f3.format(*rhi)
        rmi = numpy.linalg.norm(rlo)
        rma = numpy.linalg.norm(rhi)
        rc = numpy.linalg.norm(xc-xr)
        self.Rlimits = [rmi,rma]
        print 'findreactor.__init__ Rlimits'+f2.format(*self.Rlimits),'center-to-center distance'+f1.format(rc)

        # min,max azimuth
        phi1 = numpy.arctan2(xmi[1]-xr[1],xma[0]-xr[0])
        phi2 = numpy.arctan2(xma[1]-xr[1],xmi[0]-xr[0])
        self.philimits = [phi1,phi2]
        print 'findreactor.__init__ phi1,phi2'+f2.format(phi1,phi2)

        # min,max cos(theta)
        ctmi,ctma = 2., -2.
        for dz in [xmi[2]-xr[2],xma[2]-xr[2]]:
            for dy in [xmi[1]-xr[1],xma[1]-xr[1]]:
                for dx in [xmi[0]-xr[0],xma[0]-xr[0]]:
                    ct = dz/math.sqrt(dz*dz+dy*dy+dx*dx)
                    ctmi = min(ct,ctmi)
                    ctma = max(ct,ctma)
        self.ctlimits = [ctmi,ctma]
        print 'findreactor.__init__ ctmi,ctma'+f2.format(ctmi,ctma)
        
        return
    def main(self):

        debug = self.debug

        
        if False:  # checks R^2 distribution generation
            N = 100000
            s = self.genRm2(N=N,Rmin=self.Rlimits[0])
            s = s[s<self.Rlimits[1]]
            self.plotRm2(s,Rmin=self.Rlimits[0],Rmax=self.Rlimits[1])
            sys.exit('findreactor.main Checked R^2 distribution generation')

            
# calculate the center of every module
        centers,cX,cY,cZ = self.calcMcenters()
        self.cX,self.cY,self.cZ = cX,cY,cZ
        if debug>0:
            print 'findreactor.main centers',
            for key in sorted(centers):
                print '{0} {1:.1f} {2:.1f} {3:.1f}'.format(key,centers[key][0],centers[key][1],centers[key][2]),
            print '\n'

        showSA = True
        # calculate approximate solid angle of each module
        if showSA :
            solidA = self.calcSolidA(cX,cY,cZ,self.Xreactor)
            words = 'x {0:.1f} y {1:.1f} z {2:.1f}'.format(*self.Xreactor)
            self.plotForD(solidA,'solidAngle nominal '+words)
            for Xr in self.XreactorAlt:
                solidA = self.calcSolidA(cX,cY,cZ,Xr)
                words = 'x {0:.1f} y {1:.1f} z {2:.1f}'.format(*Xr)
                self.plotForD(solidA,'solidAngle alternative '+words)
                
            


        TotalExpt = self.TotalExpt
        Solutions = {} # {expt# : [Nibd,chi2,param]}
        N = self.GenEvtsPerExpt # starting events per expt
        
# generate the hits            
        uni = False
        incenters = False
        if uni or incenters : TotalExpt = 1

        for Nexpt in range(TotalExpt):

            self.firstTime = firstTime = Nexpt==1

            if uni:
                N = 100000
                print '\nfindreactor.main *** TESTING WITH A UNIFORM DISTRIBUTION *** \n'
            if incenters :
                print '\nfindreactor.main ++++++ TESTING WITH 1 HIT PER MODULE ++++++ \n'
                x,y,z = [],[],[]
                for k in centers:
                    q = centers[k]
                    x.append(q[0])
                    y.append(q[1])
                    z.append(q[2])
                X,Y,Z = numpy.array(x),numpy.array(y),numpy.array(z)
            else:
                X,Y,Z = self.genSP(N, uniform=uni)

    # process the hits
            X,Y,Z = self.checkSP(X,Y,Z)
            if firstTime : self.plotHits(X,Y,Z)
            iX,iY,iZ,iM = self.getModule(X,Y,Z)
            if self.debug > 0:
                self.plotHvM(iX,X,self.nx,'X')
                self.plotHvM(iY,Y,self.ny,'Y')
                self.plotHvM(iZ,Z,self.nz,'Z')
            if firstTime : self.plotModule(iX,iY,iZ,iM)

            Mfreq,Nhit = self.getMcount(iM,words='Module')
            Nibd = sum(Mfreq.values())
            if self.debug > 1 :
                self.getMcount(iX,words='iX')
                self.getMcount(iY,words='iY')
                self.getMcount(iZ,words='iZ')

            self.Nhit = Nhit

            # compare generated distribution with calculation
            checkGen = firstTime
            if checkGen : 
                self.checkU()
            
            res,method = self.findR()
            success = res['success']
            fitParam= res['x']
            chi2    = res['fun']
            Solutions[Nexpt] = [Nibd,chi2]
            Solutions[Nexpt].extend(fitParam)
            #print 'findreactor.main Solutions',Solutions

            checkStart = True
            if checkStart:
                for xstart in self.startPoints:
                    resC,methodC = self.findR(xstart=xstart,spew=False)
                    successC = res['success']
                    fitParamC= res['x']
                    chi2C    = res['fun']
                    fC = ' {0:.1f} {1:.1f} {2:.1f}'
                    words = '{0:>5} {1:>8} {2:.1f}'.format(Nexpt,methodC,successC)
                    words+= fC.format(*xstart)
                    words+= fC.format(*fitParamC[1:])
                    print 'findreactor.main checkStart',words

            scanIt = False
            if scanIt and TotalExpt <=1 : 
                if success :
                    self.scanR(plot=True, inputParam=fitParam)
                    self.scanR(plot=True, inputParam=None)
                else:
                    self.scanR(plot=True, inputParam=None)

#  summarize results
        print '{0:>5} {1:>8} chi2 C x0 y0 z0 method {2}'.format('Exp#','N(IBD)',method)
        for Nexpt in range(TotalExpt):
            f = '{0:>5} {1:>8} {2:.1f} {3:.3g} {4:.1f} {5:.1f} {6:.1f}'
            if self.nFitPar==3: f = '{0:>5} {1:>8} {2:.1f} {3:.3g} {4:.1f} {5:.1f}'
            if self.nFitPar==2: f = '{0:>5} {1:>8} {2:.1f} {3:.3g} {4:.1f}'
            print f.format(Nexpt,*Solutions[Nexpt])
        self.plotSolutions(TotalExpt,Solutions,method)
        return
    def checkU(self):
        '''
        compare Nhit distribution with calculated distribution
        '''
        Nhit = self.Nhit
        param = self.initparam(xr=self.Xreactor)
        U = self.uR2(param)
        self.plotFandD(Nhit,U,words='checkU')
        return        
    def plotSolutions(self,TotalExpt,Solutions,method):
        '''
        plots of Solutions[Expt#] = [#IBD, chi2, C, x0, y0, z0]

        some desperate stuff to reduce # of x tick marks for 1d hists
        '''
        fig = plt.figure(figsize=(10,10))


        NIBD,chi2,C,x0,y0,z0 = [],[],[],[],[],[]
        labels = ['# IBD','chi2', 'C', 'x (mm)', 'y (mm)', 'z (mm)']
        for Nexpt in range(TotalExpt):
            S = Solutions[Nexpt]
            NIBD.append( S[0] )
            chi2.append( S[1] )
            C.append(    S[2] )
            x0.append(   S[3] )
            if self.nFitPar==4:
                y0.append( S[4] )
                z0.append( S[5] )
            elif self.nFitPar==3:
                y0.append( S[4] )
                z0.append( self.Xreactor[2] )
            elif self.nFitPar==2:
                y0.append( self.Xreactor[1] )
                z0.append( self.Xreactor[2] )
                
        bigList = [NIBD, chi2, C, x0]
        nx,ny = 3,3
        if self.nFitPar==2:
            nx,ny = 2,2
        elif self.nFitPar==3:
            bigList.append( y0 )
            ny = 2
        if self.nFitPar==4:
            bigList.extend( [y0,z0] )
            nx = 3
                
        for i,q in enumerate(bigList):
            av,sd = numpy.mean(q),numpy.std(q)
            lab = ' $\mu$ {0:.1f} $\sigma$ {1:.1f}'.format(av,sd)
            plt.subplot(ny,nx,i+1) # y by x
            plt.hist(q)
            a = plt.title(labels[i] + lab,fontsize=12)
            u,v = plt.xticks()
            w = [z for z in u]
            del w[::2]
            u,v = plt.xticks(w)

        if self.nFitPar>=3:
            i += 1
            plt.subplot(ny,nx,i+1)
            plt.title('y vs x')
            h = plt.hist2d(x0,y0,cmin=0.)
            plt.colorbar(h[3])
            if self.nFitPar==4:
                i += 1
                plt.subplot(3,3,i+1)
                plt.title('y v z')
                h = plt.hist2d(z0,y0,cmin=0.)
                plt.colorbar(h[3])
                i += 1
                plt.subplot(3,3,i+1)
                plt.title('x v z')
                h = plt.hist2d(z0,x0,cmin=0.)
                plt.colorbar(h[3])

        # add global title
        words = 'Reactor x {0:.1f} y {1:.1f} z {2:.1f} '.format(*self.Xreactor)
        words += ' Nexpt {0} method {1}'.format(TotalExpt,method)
        n3 = ' {0:0.0f} {1:0.0f} {2:0.0f}'
        words += ' nxyz'+n3.format(*self.nxyz)

        plt.suptitle(words, fontsize=14)
        fig.text(0.50, 0.015,self.now,ha='center')
        fig.tight_layout(rect=[0., 0.03,1.,0.97])  # does this solve the pb of axis label overlap?

        self.showOrDraw('solutions')
        return
    def showOrPlot(self,words):
        self.showOrDraw(words)
        return
    def showOrDraw(self,words):
        '''
        show plot interactively or draw to file

        also set number of ticks 
        and clear plot after showing or drawing
        '''

        if self.drawToFile:
            filename = self.titleAsFilename(words)
            pdf = self.figDir + '/' + filename + '.pdf'
            plt.savefig(pdf)
            print 'findreactor.showOrDraw Wrote',pdf
        else:
            plt.show()
        plt.clf()
        return
    def titleAsFilename(self,title):
        '''
        return ascii suitable for use as a filename
        list of characters to be replaced is taken from https://stackoverflow.com/questions/4814040/allowed-characters-in-filename
        '''
        r = {'_': [' ', ',',  '\\', '/', ':', '"', '<', '>', '|'], 'x': ['*']}
        filename = title
        filename = ' '.join(filename.split()) # substitutes single whitespace for multiple whitespace
        for new in r:
            for old in r[new]:
                if old in filename : filename = filename.replace(old,new)
        return filename
    def scanR(self,plot=False,inputParam=None):
        '''
        scan to find reactor position
        if plot, then show some 2d plots
        '''
        if inputParam is not None:
            param0 = inputParam
            if self.nFitPar==2 : param0.extend(self.Xreactor[1:3])
            if self.nFitPar==3 : param0.append(self.Xreactor[2])
            words = 'about best fit'
        else:
            param0 = self.initparam(xr=self.Xreactor)
            words = 'about default reactor position' 


        Nhit = self.Nhit

        dz = 0.
        results = []
        x1,x2,xstep = -200.,200.,200.
        x2 += xstep
        y1,y2,ystep = -200.,200.,200.
        y2 += ystep
        print 'findreactor.scanR xstep {0:.1f} ystep {1:.1f} '.format(xstep,ystep) + words
        funmin = 1.e20
        for dx in numpy.arange(x1,x2,xstep):
            for dy in numpy.arange(y1,y2,ystep):
                xr = [param0[1]+dx,param0[2]+dy,param0[3]+dz]
                param = self.initparam(xr=xr)
                fun = self.funcR2(param)
                if plot:
                    U = self.uR2(param)
                    words = 'x {0:.1f} y {1:.1f} z {2:.1f} chi2 {3:.1f}'.format(param[1],param[2],param[3],fun)
                    self.plotFandD(Nhit,U,words)

                if self.debug > 1 : print 'findreactor.scanR fun,param',fun,param
                a = [fun]
                a.extend(param)
                results.append(a)
                funmin = min(funmin,fun)
        sr = sorted(results, key=lambda x: x[0])
        print 'Chi2 C x0 y0 z0'
        for r in sr:
            print '{0} {1:.4g} {2:.1f} {3:.1f} {4:.1f}'.format(*r)
        return
    def findR(self,xstart=None,spew=True):
        '''
        minimization to find the reactor position
        return res dict from minimizer and method

        Study with 100k IBD/expt and 10 expts shows that only Powell and SLSQP converge
        '''
        bestMethod = 'Powell' # 'SLSQP'
        studyMethods = False
        if studyMethods:
            METHODS = ['Nelder-Mead','Powell','CG','BFGS','Anneal','COBYLA','SLSQP'] # no constraints
            jacobMETHODS = ['Newton-CG'] # requires jacobian
            otherMETHODS = ['L-BFGS-B','TNC','COBYLA','SLSQP'] # allow constraints
        else:
            METHODS = [bestMethod]
            
        RES = {}
        for method in METHODS:
            param = self.initparam(xr=xstart)
            if self.firstTime :
                a = ' {:.3g}'.format(param[0])
                for x in param[1:]: a += ' {:.1f}'.format(x)
                print 'findreactor.findR initial param',a
            res = minimize(self.funcR2, param, method=method)
            RES[method] = res
            #print 'findreactor.findR method,res',method,res
            p1 = 'findreactor.findR method '+method
            for word in ['status','success','message','nit']:
                if word in res: p1 += ' {0} {1}'.format(word,res[word])
            p1 +=' chi2 {0:.1f} param'.format(res['fun'])
            for i,q in enumerate(res['x']):
                f = ' {:.1f}'
                if i==0 : f = ' {:.3g}'
                p1 += f.format(q) #            {1:.3g} {2:.1f} {3:.1f} {4:.1f}'.format(res['fun'],*res['x'])
            if spew : print p1
        res = RES[bestMethod]
        return res,bestMethod 
    def initparam(self,xr=None):
        '''
        return param = initialized parameters for funcR2
        set C = param[0] to give the count in module #0
        set param[1:] = reactor position = input
        Note that length of param depends on self.nFitPar
        '''
        cX,cY,cZ,Nhit = self.cX,self.cY,self.cZ,self.Nhit
        
        if xr is None : xr = numpy.array(self.Xreactor)

        if self.debug > 1 : print 'findreactor.initparam xr',xr

        R2 = (cX[0]-xr[0])*(cX[0]-xr[0])
        R2+= (cY[0]-xr[1])*(cY[0]-xr[1])
        R2+= (cZ[0]-xr[2])*(cZ[0]-xr[2])
        C = Nhit[0]*R2
        param = [C, xr[0]]
        if self.nFitPar==3: param.append( xr[1] )
        if self.nFitPar==4: param.extend( [xr[1],xr[2]] )
        return param
    def funcR2(self,param):
        '''
        function to be minimized. Either simple chi2 or loglike approach (eqn 39.16)

        See Review of Particle Properties (2016) eqn 39.16 

        funcR2 = 2 sum_i (u_i - n_i + n_i * log(n_i/u_i))
        where 
        u_i = C/((Xi-x0)**2 + (Yi-y0)**2 + (Zi-z0)**2) and
        n_i = # of counts
        i runs over modules
        C = param[0]
        x0,y0,z0 = param[1],param[2],param[3]
        '''
        U = self.uR2(param)
        Nhit = self.Nhit
        if self.useChi2:
            if not U.all(): print 'findreactor.funcR2 usingChi2 >1 element of U=0',U
            fun = sum((U-Nhit)*(U-Nhit)/U)
        else:
            fun = 2*sum(U - Nhit + Nhit*numpy.log(Nhit/U,where=Nhit>0.))
        return fun
    def uR2(self,param):
        '''
        uR2 = calculated rate in each module, in increasing module number
        Deal with case where number of input parameters is <4
        '''
        C,x0 = param[0],param[1]
        if self.nFitPar==4:
            y0 = param[2]
            z0 = param[3]
        elif self.nFitPar==3:
            y0 = param[2]
            z0 = self.Xreactor[2]
        elif self.nFitPar==2:
            y0 = self.Xreactor[1]
            z0 = self.Xreactor[2]
        else:
            sys.exit('findreactor.uR2 ERROR self.nFitPar '+str(self.nFitPar))
        XYZ = [x0,y0,z0]

        if self.debug>2 : print 'findreactor.uR2 param',param,'x0,y0,z0',x0,y0,z0,'XYZ',XYZ
            
        cX,cY,cZ = self.cX,self.cY,self.cZ
        R2 = (cX-x0)*(cX-x0) + (cY-y0)*(cY-y0) + (cZ-z0)*(cZ-z0)
        solidA = self.calcSolidA(cX,cY,cZ,XYZ)
        R2 *= solidA[0]/solidA
        U  = C/R2 
        return U       
    def getMcount(self,iM,words='Module'):
        '''
        get count and frequency for each module
        return dict of frequency per module and array of frequency per module ordered by increasing module number
        '''
        freq = collections.Counter(iM)
        tot = sum(freq.values())

        nhit = []
        for i in range(0,self.totalModules):
            c = 0
            if i in freq: c = freq[i]
            nhit.append(c)
        Nhit = numpy.array(nhit)
        if self.firstTime : print 'findreactor.getMcount',words,'len(iM),tot,len(freq)',len(iM),tot,len(freq)
        if self.debug>0: print 'findreactor.getMcount freq',freq
        return freq,Nhit
    def plotHvM(self,iV,V,nv,vname):
        '''
        plot module number vs position if there are more than 1 bins
        '''
        if nv==1:
            print 'findreactor.plotHvM',vname,'iV,V,nv',iV,V,nv
        else:
            h = plt.hist2d(iV,V,bins=[nv,nv],cmin=0)
            plt.title(vname)
            plt.colorbar(h[3])
            self.showOrPlot(vname)
        return
    def plotForD(self,A,words='nothing'):
        '''
        2d plot (Y v X) of fit or data i
        '''
        h = plt.hist2d(self.cX,self.cY,bins=[self.nx,self.ny],weights=A)
        plt.title(words)
        plt.colorbar(h[3])
        self.showOrPlot(words)
#        plt.show()
        return
    def plotFandD(self,D,F,words=''):
        '''
        2d plots (Y vs X) of data, fit and data-fit
        and data/fit if there are no zeros in fit
        '''
        
        torder = ['Data ','Fit ','Data-Fit ']
        t = [D,F,D-F]
        if F.all() :
            torder.append('DataOverFit ')
            t.append(D/F)

        plt.figure(figsize=(10,10))
        
        for i,q in enumerate(t):
            plt.subplot(2,2,i+1)
            h = plt.hist2d(self.cX,self.cY,bins=[self.nx,self.ny],weights=q)
            plt.title(torder[i]+words)
            plt.colorbar(h[3])
        self.showOrPlot('FandD '+words)
        #plt.show()
        return
        
    def plotHits(self,X,Y,Z):
        '''
        plot some distributions of hits
        '''
        plotR = False
        plotZvs = False
        oneD  = False
        finer = False
        if plotR : R = numpy.sqrt(X*X + Y*Y + Z*Z)
        xfine = [1]
        if finer: xfine.append(10)
        for xf in xfine:
            cf = ''
            if xf!=1 : cf = ' '+ str(xf) + 'xfiner bins'
            if oneD:
                plt.hist(X,bins=self.nx*xf,histtype='step')
                plt.title('Xhits'+cf)
                self.showOrPlot('Xhits'+cf)
#                plt.show()
                plt.hist(Y,bins=self.ny*xf,histtype='step')
                plt.title('Yhits'+cf)
                self.showOrPlot('Yhits'+cf)
#                plt.show()
            if plotR:
                plt.hist(R,20*xf,histtype='step')
                plt.title('Rhits'+cf)
                self.showOrPlot('Rhits'+cf)
                #plt.show()
                
            h = plt.hist2d(X,Y,bins=[self.nx*xf,self.ny*xf],cmin=0)
            plt.title('Yhits vs Xhits'+cf)
            plt.colorbar(h[3])
            self.showOrPlot('Yhits vs Xhits'+cf)
            #plt.show()

            if plotZvs:
                h = plt.hist2d(X,Z,bins=[self.nx*xf,self.nz*xf],cmin=0)
                plt.title('Zhits vs Xhits'+cf)
                plt.colorbar(h[3])
                self.showOrPlot('Zhits vs Xhits'+cf)
                #plt.show()

                h = plt.hist2d(Y,Z,bins=[self.ny*xf,self.nz*xf],cmin=0)
                plt.title('Zhits vs Yhits'+cf)
                plt.colorbar(h[3])
                self.showOrPlot('Zhits vs Xhits'+cf)
                #plt.show()
        
        return
    def plotModule(self,jX,jY,jZ,iM):
        '''
        plot iX, iY, iY vx iX and maybe other stuff
        '''
        plot1D = False
        delta = 0
        iX,iY,iZ = jX+delta,jY+delta,jZ+delta
        if plot1D:
            plt.hist(iX,bins=self.nx,histtype='step')
            plt.title('Xmodule')
            self.showOrPlot('Xmodule')
            #plt.show()
            plt.hist(iY,bins=self.ny,histtype='step')
            plt.title('Ymodule')
            self.showOrPlot('Ymodule')
            #plt.show()
        h = plt.hist2d(iX,iY,bins=[self.nx,self.ny],cmin=0)
        if self.debug > 2: print 'findreactor.plotModule hist2d returns h[0]',h[0]
        plt.title('Ymodule vs Xmodule')
        plt.colorbar(h[3])
        self.showOrPlot('Ymodule vs Xmodule')
        #plt.show()
        return
    def getModNum(self,ix,iy,iz):
        ''' define module number in one location '''
        return ix + iy*self.nx + iz*self.ny*self.nx
    def getModule(self,X,Y,Z):
        '''
        return module numbers and global module number given spacepoints
        '''
        ix = numpy.floor((X - self.minXdet[0])/self.dXmod[0]).astype(int)
        iy = numpy.floor((Y - self.minXdet[1])/self.dXmod[1]).astype(int)
        iz = numpy.floor((Z - self.minXdet[2])/self.dXmod[2]).astype(int)
        #im = iz*self.ny*self.nx + iy*self.nx + ix # global module number
        im = self.getModNum(ix,iy,iz)
        xfail = min(ix)<0 or max(ix)>=self.nx 
        yfail = min(iy)<0 or max(iy)>=self.ny 
        zfail = min(iz)<0 or max(iz)>=self.nz
        mfail = min(im)<0 or max(im)>=self.nz*self.ny*self.nx
        if xfail or yfail or zfail or mfail : 
            print 'findreactor.getModule xfail,yfail,zfail,mfail', xfail,yfail,zfail,mfail
            print 'findreactor.getModule min(ix),max(ix),min(iy),max(iy),min(iz),max(iz),min(im),max(im)',min(ix),max(ix),min(iy),max(iy),min(iz),max(iz),min(im),max(im)
            
        return ix,iy,iz,im
    def plotRm2(self,s,Rmin=0.,Rmax=None):
        '''
        plot input distribution s in ordinate range [Rmin,Rmax]
        print number of entries <Rmin, >Rmax, and in range [Rmin,Rmax]
        '''
        if Rmax is None or Rmax<Rmin:
            counts, bins, _ = plt.hist(s, bins=20, histtype='step')
        else:
            counts, bins, _ = plt.hist(s, range=(Rmin,Rmax), bins=20, histtype='step' )
        lo = sum(s<Rmin)
        hi = sum(s>Rmax)
        inrange = len(s)-lo-hi
        print 'findreactor.plotRm2 lo,inrange,hi',lo,inrange,hi
        plt.title('supposed to follow 1/x^2')
        self.plotOrShow('supposed to follow 1/x^2')
#        plt.show()
        return
    def genSP(self,N=1,uniform=False):
        '''
        generate up to N space points 
        if uniform is True, just uniformly distribute in x,y,z
        else try to generate according to 1/R2 nearly inside the detector
        '''
        if uniform:
            X  = self.genRange(N,V=[self.minXdet[0],self.maxXdet[0]])
            Y  = self.genRange(N,V=[self.minXdet[1],self.maxXdet[1]])
            Z  = self.genRange(N,V=[self.minXdet[2],self.maxXdet[2]])
        else:
            R = self.genRm2(N,Rmin=self.Rlimits[0])
            R = R[R<self.Rlimits[1]]
            n = len(R)
            Phi = self.genRange(n,V=self.philimits)
            CT  = self.genRange(n,V=self.ctlimits)
            ST  = numpy.sin( numpy.arccos( CT ) )
            X = self.Xreactor[0] + R*ST*numpy.cos(Phi)
            Y = self.Xreactor[1] + R*ST*numpy.sin(Phi)
            Z = self.Xreactor[2] + R*CT
        return X,Y,Z
    def checkSP(self,X,Y,Z):
        '''
        count generated space points are inside the detector
        return only space points generated inside detector
        '''
        N = len(X)
        allcheck = self.inDet(X,Y,Z)
        check = sum( allcheck )
        if self.firstTime : print 'findreactor.checkSP #points',N,'#in det.',check

        X = X[allcheck]
        Y = Y[allcheck]
        Z = Z[allcheck]
        allcheck = self.inDet(X,Y,Z)
        check = sum(allcheck)
        if self.firstTime : print 'findreactor.checkSP after check. #in det.',check

        return X,Y,Z
    def inDet(self,X,Y,Z):
        '''
        return boolean array that is true if X,Y,Z combination is inside detector
        '''
        x1,x2 = self.minXdet[0],self.maxXdet[0]
        y1,y2 = self.minXdet[1],self.maxXdet[1]
        z1,z2 = self.minXdet[2],self.maxXdet[2]
        Xcheck = (x1<X)*(X<x2)
        Ycheck = (y1<Y)*(Y<y2)
        Zcheck = (z1<Z)*(Z<z2)
        allcheck = Xcheck * Ycheck * Zcheck
        return allcheck
    def calcMcenters(self):
        '''
        calculate the centers of all modules
        return centers as dict[module# : [x,y,z]] and arrays cX,cY,cZ in increasing module number order
        '''
        vmi,vma = numpy.ndarray.copy(self.minXdet),numpy.ndarray.copy(self.maxXdet)
        nxyz = numpy.ndarray.copy(self.nxyz).astype(int)
        dv = (vma-vmi)/nxyz.astype(float)
        print 'findreactor.calcMcenters vmi,vma,nxyz,dv',vmi,vma,nxyz,dv
        centers = {}
        for jz in range(0,nxyz[2]):
            z = vmi[2] + (float(jz)+0.5)*dv[2]
            for jy in range(0,nxyz[1]):
                y = vmi[1] + (float(jy)+0.5)*dv[1]
                for jx in range(0,nxyz[0]):
                    x = vmi[0] + (float(jx)+0.5)*dv[0]
#                    jm = jx + nxyz[0]*jy + nxyz[0]*nxyz[1]*jz
                    jm = self.getModNum(jx,jy,jz)
                    if jm in centers:
                        print 'findreactor.calcMcenters ERROR jm',jm,'in centers',centers[jm],'current x,y,z',x,y,z,'current jx,jy,jz',jx,jy,jz
                        sys.exit('findreactor.calcMcenters ERROR jm '+str(jm))
                    centers[jm] = [x,y,z]
                    
        keylist = sorted(centers)
        x,y,z = [],[],[]
        for k in keylist:
            q = centers[k]
            x.append(q[0])
            y.append(q[1])
            z.append(q[2])
        cX,cY,cZ = numpy.array(x),numpy.array(y),numpy.array(z)
        return centers,cX,cY,cZ
    def calcSolidA(self,cX,cY,cZ,xr):
        '''
        given centers of each module in cX,cY,cZ, and reactor position xr
        calculate approximate solid angle of all modules.
        approximation is the area of the plane perpendicular to the 
        reactor center-module center within the module divided by 4pi*R^2

        20210505 Just use the cross-section area / (4*pi*R^2) where 
        area = length * width of module
        R = distance from module center to reactor center
        '''
        if self.debug > 1 : print 'findreactor.calcSolidA module#,solidAngle'
        w,l = self.dXmod[0], self.dXmod[2]
        rx,ry,rz = xr[0],xr[1],xr[2]
        solidA = []
        for i,yy in enumerate(cX):
            dx2,dy2,dz2 = cX[i]-rx,cY[i]-ry,cZ[i]-rz
            dx2,dy2,dz2 = dx2*dx2, dy2*dy2, dz2*dz2
            sA = w*l / (dx2 + dy2 + dz2)/4./math.pi
            solidA.append(sA)
            if self.debug > 1 : 
                if i%20==0 : print '\n'
                print '{0} {1:.5f}'.format(i,sA),
        if self.debug > 1 : print '\n'
        if self.debug > 2 :
            print 'findreactor.calcSolidA module#,solidAngle,cX,cY,cZ'
            for i,sa in enumerate(solidA):
                print '{0:3} {1:.5f} {2:.1f} {3:.1f} {4:.1f}'.format(i,sa,cX[i],cY[i],cZ[i])

        return numpy.array(solidA)
    def OLDcalcSolidA(self,cX,cY,cZ,xr):
        '''
        DEPRECATED
        given centers of each module in cX,cY,cZ, and reactor position xr
        calculate approximate solid angle of all modules.
        approximation is the solid angle subtended by 
        the diagonal width from the top front corners to the 
        back bottom corners of each module
        '''
        if self.debug > 1 : print 'findreactor.OLDcalcSolidA module#,solidAngle'
        dx,dy,dz = self.dXmod
        rx,ry,rz = xr[0],xr[1],xr[2]
        solidA = []
        for i,yy in enumerate(cX):
            x,y,z = cX[i]-rx,cY[i]-ry,cZ[i]-rz
            x1,x2 = x-dx/2.,x+dx/2.
            y1,y2 = y-dy/2.,y+dz/2.
            z1,z2 = z-dz/2.,z+dz/2.
            phi1,phi2 = numpy.arctan2(y2,x1),numpy.arctan2(y1,x2)
            ct1,ct2 = 2.,-2.
            for u in [x1,x2]:
                for v in [y1,y2]:
                    for w in [z1,z2]:
                        ct = w/math.sqrt(v*v+u*u+w*w)
                        ct1=min(ct1,ct)
                        ct2=max(ct2,ct)
            sA = abs(phi1-phi2)*abs(ct1-ct2)
            solidA.append(sA)
            if self.debug > 1 : 
                if i%20==0 : print '\n'
                print '{0} {1:.5f}'.format(i,sA),
        if self.debug > 1 : print '\n'
        return numpy.array(solidA)
    def genRm2(self,N=1,Rmin=0.):
        ''' generate N samples from 1/R^2 distribution starting at Rmin'''
        a, m = 1., Rmin
        s = (numpy.random.pareto(a, N) + 1.)*m
        return s
    def genRange(self,N=1,V=[0.,1]):
        '''
        return random floats in interval [V(0),V(1))
        '''
        Vmax = max(V)
        Vmin = min(V)
        s = numpy.random.random(N) * (Vmax-Vmin) + Vmin
        return s
    def testPareto(self):
        a, m = 1., 6700. 
        x1,x2 = m,m+(9200.-m)
        s = (numpy.random.pareto(a, 1000000) + 1)*m
        count, bins, _ = plt.hist(s, bins=20, range=(x1,x2), histtype='step')#100, color='blue') #, density=True)
        #for y,x in zip(count,bins): print y,x
        dx =  bins[1]-bins[0]
        bins = bins + dx/2.
        fit = a*m**a / bins**(a+1)
        plt.plot(bins, max(count)*fit/max(fit), linewidth=2, color='r')

        fit2 = count[0]*(bins[0]/bins)**2
        plt.plot(bins, fit2, linewidth=2, color='black')
        
        
        #plt.yscale('log')
        self.plotOrShow('testPareto')
#        plt.show()
        return
if __name__ == '__main__' :
    debug = -1
    TotalExpt = 10
    drawToFile = False
    nFitPar = 4
    GenEvtsPerExpt = 1000000
    useChi2 = False
    setNZ = 1
    if len(sys.argv)>1 :
        w = sys.argv[1]
        if 'help' in w.lower():
            print 'USAGE: python findreactor.py [debug] [TotalExpt] [drawToFile] [nFitPar] [GenEvtsPerExpt] [useChi2] [setNZ]'
            print 'DEFAULTS: python findreactor.py',debug,TotalExpt,drawToFile,nFitPar,GenEvtsPerExpt,useChi2,setNZ
            sys.exit('help was provided')
    if len(sys.argv)>1 : debug = int(sys.argv[1])
    if len(sys.argv)>2 : TotalExpt = int(sys.argv[2])
    if len(sys.argv)>3 : drawToFile = bool(sys.argv[3])
    if len(sys.argv)>4 : nFitPar = int(sys.argv[4])
    if len(sys.argv)>5 : GenEvtsPerExpt = int(sys.argv[5])
    if len(sys.argv)>6 : useChi2 = bool(sys.argv[6])
    if len(sys.argv)>6 : setNZ = max(1,int(sys.argv[7]))
    nFitPar = max(2,min(4,nFitPar))
    fr = findreactor(debug,TotalExpt=TotalExpt,drawToFile=drawToFile,nFitPar=nFitPar,GenEvtsPerExpt=GenEvtsPerExpt,useChi2=useChi2,setNZ=setNZ)
    fr.main()

                
            
    
